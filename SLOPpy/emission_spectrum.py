from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ['compute_emission_spectrum',
           'plot_emission_spectrum',
           'compute_emission_spectrum_planetRF',
           'plot_emission_spectrum_planetRF',
           'compute_emission_spectrum_stellarRF',
           'plot_emission_spectrum_stellarRF',
           'compute_emission_spectrum_observerRF',
           'plot_emission_spectrum_observerRF',
           'compute_emission_spectrum_core',
           'plot_emission_spectrum_core',
           'weighted_average_spectra',
           'emission_template_from_range',
           'average_and_bin',
           'emission_ylim']


subroutine_name = 'emission_spectrum'


def compute_emission_spectrum_planetRF(config_in, lines_label):
    compute_emission_spectrum(config_in, lines_label, reference='planetRF')


def plot_emission_spectrum_planetRF(config_in, lines_label, night_input=''):
    plot_emission_spectrum(config_in, lines_label, night_input, reference='planetRF')


def compute_emission_spectrum_stellarRF(config_in, lines_label):
    compute_emission_spectrum(config_in, lines_label, reference='stellarRF')


def plot_emission_spectrum_stellarRF(config_in, lines_label, night_input=''):
    plot_emission_spectrum(config_in, lines_label, night_input, reference='stellarRF')


def compute_emission_spectrum_observerRF(config_in, lines_label):
    compute_emission_spectrum(config_in, lines_label, reference='observerRF')


def plot_emission_spectrum_observerRF(config_in, lines_label, night_input=''):
    plot_emission_spectrum(config_in, lines_label, night_input, reference='observerRF')


def compute_emission_spectrum(config_in, lines_label, reference='planetRF'):

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    compute_emission_spectrum_core(config_in, subroutine_name, lines_dict['range'],
                                   reference=reference, lines_label=lines_label)


def plot_emission_spectrum(config_in, lines_label, night_input='', reference='planetRF'):

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    plot_range = lines_dict.get('plot_range', lines_dict['range'])

    plot_emission_spectrum_core(config_in, subroutine_name, plot_range,
                                reference=reference, lines_label=lines_label, night_input=night_input)


def weighted_average_spectra(spectra, spectra_err):
    """
    Weighted average of a set of spectra, with weights equal to the inverse of the variance.
    Pixels with null or negative errors (i.e., not covered by the observations) are excluded
    :param spectra: 2D array (n_spectra x n_pixels)
    :param spectra_err: 2D array with the errors
    :return: average, average_err, sum of the weights. Pixels without data have
        average and error equal to zero
    """
    spectra = np.atleast_2d(spectra)
    spectra_err = np.atleast_2d(spectra_err)

    valid = (spectra_err > 0.) & np.isfinite(spectra) & np.isfinite(spectra_err)
    weights = np.zeros_like(spectra_err)
    weights[valid] = 1. / spectra_err[valid]**2

    sum_weights = np.sum(weights, axis=0)
    covered = (sum_weights > 0.)

    average = np.zeros(np.shape(spectra)[1], dtype=np.double)
    average_err = np.zeros(np.shape(spectra)[1], dtype=np.double)

    average[covered] = np.sum(weights[:, covered] * np.where(valid, spectra, 0.)[:, covered], axis=0) \
        / sum_weights[covered]
    average_err[covered] = 1. / np.sqrt(sum_weights[covered])

    return average, average_err, sum_weights


def emission_template_from_range(shared_data, subroutine, wave_range=None):
    """
    Wavelength grid of the emission spectrum, selected from the shared coadded and binned grids
    :param wave_range: [wave_min, wave_max]; the full shared grid is used if None
    """
    if wave_range is None:
        return {
            'subroutine': subroutine,
            'range': shared_data['coadd']['wavelength_range'],
            'wave': shared_data['coadd']['wave'],
            'step': shared_data['coadd']['step'],
            'size': shared_data['coadd']['size'],
            'binned_wave': shared_data['binned']['wave'],
            'binned_step': shared_data['binned']['step'],
            'binned_size': shared_data['binned']['size']
        }

    shared_selection = (shared_data['coadd']['wave'] >= wave_range[0]) \
        & (shared_data['coadd']['wave'] < wave_range[1])
    binned_selection = (shared_data['binned']['wave'] >= wave_range[0]) \
        & (shared_data['binned']['wave'] < wave_range[1])

    return {
        'subroutine': subroutine,
        'range': wave_range,
        'wave': shared_data['coadd']['wave'][shared_selection],
        'step': shared_data['coadd']['step'][shared_selection],
        'size': int(np.sum(shared_selection)),
        'binned_wave': shared_data['binned']['wave'][binned_selection],
        'binned_step': shared_data['binned']['step'][binned_selection],
        'binned_size': int(np.sum(binned_selection))
    }


def average_and_bin(emission, prefix, spectra, spectra_err):
    """
    Store the weighted average of the spectra and its binned version in the emission dictionary,
    under the keys 'average<prefix>', 'average<prefix>_err', 'binned<prefix>', 'binned<prefix>_err'
    """
    emission['average' + prefix], emission['average' + prefix + '_err'], emission['sum_weights' + prefix] = \
        weighted_average_spectra(spectra, spectra_err)

    emission['binned' + prefix] = \
        rebin_1d_to_1d(emission['wave'],
                       emission['step'],
                       emission['average' + prefix],
                       emission['binned_wave'],
                       emission['binned_step'],
                       preserve_flux=False)
    emission['binned' + prefix + '_err'] = \
        rebin_1d_to_1d(emission['wave'],
                       emission['step'],
                       emission['average' + prefix + '_err'],
                       emission['binned_wave'],
                       emission['binned_step'],
                       preserve_flux=False,
                       is_error=True)


def compute_emission_spectrum_core(config_in, output_name, wave_range, reference='planetRF', lines_label=''):
    """
    Emission spectrum of each eclipse night, in the selected reference frame.

    The emission spectra computed by emission_spectrum_preparation are already continuum-normalized
    and free from the blaze function, so no blaze correction, CLV+RM correction or continuum
    normalization are applied here: the spectra are only shifted to the requested reference frame,
    rebinned to 1D and averaged
    :param output_name: name of the output files
    :param wave_range: [wave_min, wave_max], or None to use the full spectrum
    """

    night_dict = from_config_get_nights(config_in, phase='eclipse')

    shared_data = load_from_cpickle('shared', config_in['output'])
    emission_template = emission_template_from_range(shared_data, output_name, wave_range)

    for night in night_dict:

        print()
        print("Running {0:45s} for  {1:20s}   Night:{2:15s}  ".format(output_name, lines_label, night))

        preparation = load_from_cpickle('emission_preparation', config_in['output'], night)

        try:
            emission = load_from_cpickle(output_name + '_' + reference, config_in['output'], night, lines_label)
            """ recompute if the preparation step has been computed again with different settings """
            if emission.get('savgol_pams', None) != preparation.get('savgol_pams', None):
                raise ValueError('Stored emission spectrum computed from a different preparation')
            print("{0:45s} Night:{1:15s}   {2:s}   {3:s}".format(output_name, night, lines_label, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}   {3:s}".format(output_name, night, lines_label, 'Computing'))

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        emission = emission_template.copy()
        emission['reference'] = reference
        emission['savgol_pams'] = preparation.get('savgol_pams', None)

        for obs in lists['observations']:

            emission[obs] = {
                'BJD': input_data[obs]['BJD'],
                'AIRMASS': input_data[obs]['AIRMASS']
            }

            if reference in ['observer', 'observerRF', 'ORF']:
                rv_shift = 0.000
            elif reference in ['stellar', 'stellarRF', 'SRF']:
                rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
            else:
                rv_shift = observational_pams[obs]['rv_shift_ORF2PRF']

            """ Rebin the 2D emission spectra to 1D. The spectra are already free from the
                blaze function, and they are dimensionless, so the flux is not preserved """
            emission[obs]['rebinned'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               preparation[obs]['deblazed'],
                               None,
                               emission['wave'],
                               emission['step'],
                               preserve_flux=False,
                               skip_blaze_correction=True,
                               rv_shift=rv_shift)

            emission[obs]['rebinned_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               preparation[obs]['deblazed_err'],
                               None,
                               emission['wave'],
                               emission['step'],
                               preserve_flux=False,
                               skip_blaze_correction=True,
                               rv_shift=rv_shift,
                               is_error=True)

            ### Small border bugfix
            if emission[obs]['rebinned_err'][0] == 0:
                emission[obs]['rebinned'][0] = emission[obs]['rebinned'][1]
                emission[obs]['rebinned_err'][0] = emission[obs]['rebinned_err'][1]

            if emission[obs]['rebinned_err'][-1] == 0:
                emission[obs]['rebinned'][-1] = emission[obs]['rebinned'][-2]
                emission[obs]['rebinned_err'][-1] = emission[obs]['rebinned_err'][-2]

        """ Average of the spectra in the full_transit list (planet signal) and in the
            out_transit list (reference spectra used for the master-out) """
        for prefix, list_name in [['', 'transit_full'], ['_out', 'transit_out']]:
            average_and_bin(emission, prefix,
                            [emission[obs]['rebinned'] for obs in lists[list_name]],
                            [emission[obs]['rebinned_err'] for obs in lists[list_name]])

        save_to_cpickle(output_name + '_' + reference, emission, config_in['output'], night, lines_label)

        # Forcing memory deallocation
        emission = None


def emission_ylim(values, errors=None):
    """ symmetric limits around zero, including most of the data points """
    values = np.asarray(values)
    selection = np.isfinite(values) & (values != 0.)
    if errors is not None:
        selection &= (np.asarray(errors) > 0.)
    if np.sum(selection) == 0:
        return -0.01, 0.01
    limit = 1.2 * np.percentile(np.abs(values[selection]), 99.)
    return -limit, limit


def plot_emission_spectrum_core(config_in, output_name, plot_range, reference='planetRF', lines_label='', night_input=''):

    night_dict = from_config_get_nights(config_in, phase='eclipse')

    if night_input == '':
        night_list = night_dict
    else:
        """ nights with phase: transit are skipped """
        night_list = [night for night in np.atleast_1d(night_input) if night in night_dict]

    os.system('mkdir -p plots')

    interactive_plots = from_config_get_interactive_plots(config_in)

    for night in night_list:

        filename_rad = output_name + '_' + reference

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the analysis"""
        try:
            emission = load_from_cpickle(filename_rad, config_in['output'], night, lines_label)
        except (FileNotFoundError, IOError):
            print()
            print("No emission spectrum in {0:s}, no plots".format(reference))
            continue

        """ Creation of the color array, based on the BJD of the observations
        """
        bjd = [emission[obs]['BJD'] - 2450000.0 for obs in lists['observations']]

        color_cmap = plt.cm.viridis
        color_norm = plt.Normalize(vmin=bjd[0], vmax=bjd[-1])

        fig = plt.figure(figsize=(12, 6))

        gs = GridSpec(2, 2, width_ratios=[50, 1])
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for ax, list_name in [[ax1, 'transit_full'], [ax2, 'transit_out']]:
            for obs in lists[list_name]:
                color = [color_cmap(color_norm(emission[obs]['BJD'] - 2450000.0))[:-1]]
                ax.scatter(emission['wave'],
                           emission[obs]['rebinned'],
                           c=color, s=1, zorder=3, alpha=0.25)

        ax1.set_ylim(emission_ylim(np.concatenate([emission[obs]['rebinned'] for obs in lists['transit_full']])))
        ax1.set_xlim(plot_range[0], plot_range[1])
        ax2.set_xlabel(r'$\lambda$ [$\AA$]')
        ax1.set_title('Lines: {0:s} Night: {1:s} \n Emission spectra (full_transit list) in {2:s}'.format(
            lines_label, night, reference))
        ax2.set_title('Emission spectra (out_transit list) in {0:s}'.format(reference))

        sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)

        output_file = get_filename(filename_rad + '_observations',
                                   config_in['output'], night, lines_label, extension='.pdf')
        plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
        if interactive_plots:
            plt.show()
        plt.close()

        fig = plt.figure(figsize=(12, 6))

        gs = GridSpec(2, 1)
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)

        for ax, prefix, label in [[ax1, '', 'full_transit'], [ax2, '_out', 'out_transit']]:
            ax.errorbar(emission['wave'],
                        emission['average' + prefix],
                        yerr=emission['average' + prefix + '_err'],
                        fmt='ko', ms=1, zorder=5, alpha=0.25, label='average')
            ax.errorbar(emission['binned_wave'],
                        emission['binned' + prefix],
                        yerr=emission['binned' + prefix + '_err'],
                        fmt='ro', ms=4, lw=2, zorder=10, label='binned average')
            ax.axhline(0.0, c='C0', zorder=1)

        ax1.set_ylim(emission_ylim(emission['average'], emission['average_err']))
        ax1.set_xlim(plot_range[0], plot_range[1])
        ax2.set_xlabel(r'$\lambda$ [$\AA$]')
        ax2.legend(loc=3)
        ax1.set_title('Lines: {0:s} Night: {1:s} \n Average emission spectrum (full_transit list) in {2:s}'.format(
            lines_label, night, reference))
        ax2.set_title('Average emission spectrum (out_transit list) in {0:s}'.format(reference))
        fig.subplots_adjust(hspace=0.4)

        output_file = get_filename(filename_rad + '_binned',
                                   config_in['output'], night, lines_label, extension='.pdf')
        plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
        if interactive_plots:
            plt.show()
        plt.close()
