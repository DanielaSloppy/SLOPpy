from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = [
           'compute_transmission_spectrum_average',
           'plot_transmission_spectrum_average'
           ]


subroutine_name = 'transmission_spectrum_average'
pick_files =  'transmission_spectrum'
sampler = 'emcee'


def compute_transmission_spectrum_average(config_in, lines_label, reference='planetRF', pca_iteration=-1):

    night_dict = from_config_get_nights(config_in, phase='transit')

    spectral_lines = from_config_get_spectral_lines(config_in)
    line_iter_dict = spectral_lines[lines_label]

    shared_data = load_from_cpickle('shared', config_in['output'])

    total_n_transit_full = 0
    total_n_transit_out = 0

    total_lists = {}

    """ Using the user defined range to define the transmission spectrum region 
        This range can be larger than the one defined for the MCMC range, and it
        MUST include the continuum windws for the transmission lightcurve
    """

    shared_selection = (shared_data['coadd']['wave'] >= line_iter_dict['range'][0]) \
        & (shared_data['coadd']['wave'] < line_iter_dict['range'][1])
    binned_selection = (shared_data['binned']['wave'] >= line_iter_dict['range'][0]) \
        & (shared_data['binned']['wave'] < line_iter_dict['range'][1])

    transmission_average_template = {
        'subroutine': subroutine_name + '_' + reference,
        'range': line_iter_dict['range'],
        'wave': shared_data['coadd']['wave'][shared_selection],
        'step': shared_data['coadd']['step'][shared_selection],
        'size': int(np.sum(shared_selection)),
        'binned_wave': shared_data['binned']['wave'][binned_selection],
        'binned_step': shared_data['binned']['step'][binned_selection],
        'binned_size': int(np.sum(binned_selection))
    }

    results_list = ['user',
                    'mcmc_night_MED',
                    'mcmc_night_MAP',
                    'mcmc_global_MED',
                    'mcmc_global_MAP']

    for results_selection in results_list:

        """ First check to see if we need to compute the average transmission
        iteratively when PCA has been employed """

        for night in night_dict:
            preparation_input = load_from_cpickle('transmission_preparation', config_in['output'], night)
            if preparation_input.get('pca_output', False):
                if pca_iteration >= 0:
                    it_string = str(pca_iteration).zfill(2)
                else:
                    it_string = str(pca_parameters.get('ref_iteration')).zfill(2)
            else:
                it_string = ''
            preparation = None
            break


        transmission_average = transmission_average_template.copy()

        try:
            transmission_average = load_from_cpickle(subroutine_name + '_' + reference + '_' + results_selection, config_in['output'], lines=lines_label, it_string=it_string)
            print("{0:45s}    {1:s}   {2:s}".format(subroutine_name + '_' + reference, results_selection, 'Retrieved'))
            continue
        except (FileNotFoundError, IOError):
            skip_iteration = False
            #print("{0:45s}    {1:s}   {2:s}".format(subroutine_name + '_' + reference, results_selection, 'Computing'))


        #skip_iteration = False
        for night in night_dict:
            # """ Retrieving the list of observations"""

            total_lists[night] = load_from_cpickle('lists', config_in['output'], night)
            #try:
            #    transmission_average[night] = load_from_cpickle('transmission_second_telluric_'+reference, config_in['output'], night)
            #    print("   Using transmission spectra with second telluric correction for Night: {0:s}".format(night))
            #except:
            #    transmission_average[night] = load_from_cpickle('transmission_'+reference, config_in['output'], night)

            try:
                transmission_average[night] = load_from_cpickle(pick_files + '_' + reference + '_' + results_selection, config_in['output'], night, lines_label, it_string)
            except:
                skip_iteration = True
                print(skip_iteration)
            total_n_transit_full += len(total_lists[night]['transit_full'])
            total_n_transit_out += len(total_lists[night]['transit_out'])

        if skip_iteration: continue
        print("{0:45s}    {1:s}   {2:s}".format(subroutine_name + '_' + reference, results_selection, 'Computing'))

        array_average_in = np.zeros([total_n_transit_full, transmission_average['size']])
        weights_average_in = np.zeros([total_n_transit_full, transmission_average['size']])
        clvrm_average_in = np.zeros([total_n_transit_full, transmission_average['size']])
        uncorr_average_in = np.zeros([total_n_transit_full, transmission_average['size']])
        i_total_in = 0

        array_average_out = np.zeros([total_n_transit_out, transmission_average['size']])
        weights_average_out = np.zeros([total_n_transit_out, transmission_average['size']])
        i_total_out = 0

        for night in night_dict:

            for obs in total_lists[night]['transit_full']:
                array_average_in[i_total_in, :] = transmission_average[night][obs]['normalized'][:]
                weights_average_in[i_total_in, :] = 1./(transmission_average[night][obs]['normalized_err']**2.)
                clvrm_average_in[i_total_in, :] = transmission_average[night][obs]['clv_model_rebinned'][:]
                uncorr_average_in[i_total_in, :] = transmission_average[night][obs]['normalized_uncorrected'][:]
                i_total_in += 1

            for obs in total_lists[night]['transit_out']:
                array_average_out[i_total_out, :] = transmission_average[night][obs]['normalized'][:]
                weights_average_out[i_total_out, :] = 1. / (transmission_average[night][obs]['normalized_err'] ** 2.)
                i_total_out += 1

        transmission_average['average'], transmission_average['sum_weights'] = np.average(
            array_average_in, axis = 0, weights = weights_average_in, returned = True)
        transmission_average['average_err'] = 1./np.sqrt(transmission_average['sum_weights'])

        transmission_average['average_clv_model'], _ = np.average(
            clvrm_average_in, axis = 0, weights = weights_average_in, returned = True)

        transmission_average['average_uncorrected'], _ = np.average(
            uncorr_average_in, axis = 0, weights = weights_average_in, returned = True)

        transmission_average['binned'] = \
            rebin_1d_to_1d(transmission_average['wave'],
                        transmission_average['step'],
                        transmission_average['average'],
                        transmission_average['binned_wave'],
                        transmission_average['binned_step'],
                        preserve_flux=False)
        transmission_average['binned_err'] = \
            rebin_1d_to_1d(transmission_average['wave'],
                        transmission_average['step'],
                        transmission_average['average_err'],
                        transmission_average['binned_wave'],
                        transmission_average['binned_step'],
                        preserve_flux=False,
                        is_error=True)

        transmission_average['binned_clv_model'] = \
            rebin_1d_to_1d(transmission_average['wave'],
                        transmission_average['step'],
                        transmission_average['average_clv_model'],
                        transmission_average['binned_wave'],
                        transmission_average['binned_step'],
                        preserve_flux=False)

        transmission_average['binned_uncorrected'] = \
            rebin_1d_to_1d(transmission_average['wave'],
                        transmission_average['step'],
                        transmission_average['average_uncorrected'],
                        transmission_average['binned_wave'],
                        transmission_average['binned_step'],
                        preserve_flux=False)

        transmission_average['average_out'], transmission_average['sum_weights_out'] = np.average(
            array_average_out, axis=0, weights=weights_average_out, returned=True)
        transmission_average['average_out_err'] = 1./np.sqrt(transmission_average['sum_weights_out'])

        transmission_average['binned_out'] = \
            rebin_1d_to_1d(transmission_average['wave'],
                        transmission_average['step'],
                        transmission_average['average_out'],
                        transmission_average['binned_wave'],
                        transmission_average['binned_step'],
                        preserve_flux=False)
        transmission_average['binned_out_err'] = \
            rebin_1d_to_1d(transmission_average['wave'],
                        transmission_average['step'],
                        transmission_average['average_out_err'],
                        transmission_average['binned_wave'],
                        transmission_average['binned_step'],
                        preserve_flux=False,
                        is_error=True)

        save_to_cpickle(subroutine_name + '_' + reference + '_' + results_selection, transmission_average, config_in['output'], lines=lines_label, it_string=it_string)


def plot_transmission_spectrum_average(config_in, lines_label, night_input='', results_input='', reference='planetRF', pca_iteration=-1):

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    if results_input=='':
        results_list = ['user',
                'mcmc_night_MED',
                'mcmc_night_MAP',
                'mcmc_global_MED',
                'mcmc_global_MAP']
    else:
        results_list = np.atleast_1d(results_input)

    os.system('mkdir -p plots')

    interactive_plots = from_config_get_interactive_plots(config_in)



    # Workaround to check if the transmission spectrum has been obtained through PCA iterations
    night_dict = from_config_get_nights(config_in, phase='transit')
    for night in night_dict:
        preparation_input = load_from_cpickle('transmission_preparation', config_in['output'], night)

        if preparation_input.get('pca_output', False):
            if pca_iteration >= 0:
                it_string = str(pca_iteration).zfill(2)
            else:
                it_string = str(preparation_input.get('ref_iteration')).zfill(2)
        else:
            it_string = ''
        preparation_input = None
        break



    for results_selection in results_list:

        try:
            transmission_average = load_from_cpickle(subroutine_name + '_' + reference + '_' + results_selection, config_in['output'], lines=lines_label, it_string=it_string)
            print("{0:45s}    {1:s}   {2:s}".format(subroutine_name + '_' + reference, results_selection, 'Plotting'))
        except (FileNotFoundError, IOError):
            print("{0:45s}                         {1:s}".format(subroutine_name + '_' + reference, 'Plot skipped'))
            return

        filename_rad = subroutine_name + '_' + reference + '_' + results_selection

        fig = plt.figure(figsize=(12, 9))

        gs = GridSpec(2, 1)
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)

        spec_offset = 0.025

        ax1.errorbar(transmission_average['wave'],
                    transmission_average['average'],
                    yerr=transmission_average['average_err'],
                    fmt='ko', ms=1, zorder=10, alpha=0.10, label='average')

        #ax1.scatter(transmission_average['wave'],
        #            transmission_average['average'],
        #            c='black',
        #            s=2, zorder=15,
        #            label='average',
        #            )

        ax1.errorbar(transmission_average['binned_wave'],
                    transmission_average['binned'],
                    yerr=transmission_average['binned_err'],
                    fmt='ko', ms=3, zorder=20, label='binned')

        #ax1.scatter(transmission_average['wave'],
        #            transmission_average['average'],
        #            c='black',
        #            s=2, zorder=200,
        #            label='average',
        #            )

        ax2.errorbar(transmission_average['binned_wave'],
                    transmission_average['binned_out'],
                    yerr=transmission_average['binned_out_err'],
                    fmt='ko', ms=3, zorder=20, label='binned out')

        ax2.errorbar(transmission_average['wave'],
                    transmission_average['average_out'],
                    yerr=transmission_average['average_out_err'],
                    fmt='ko', ms=1, zorder=10, alpha=0.10, label='average out')

        for n_night, night in enumerate(night_dict):
            ax1.errorbar(transmission_average['wave'],
                        transmission_average[night]['average']-spec_offset*(1.+n_night),
                        yerr=transmission_average[night]['average_err'],
                        color='C'+repr(n_night),
                        fmt='o', ms=1, zorder=1, alpha=0.25)

            ax1.scatter(transmission_average['wave'],
                        transmission_average[night]['average']-spec_offset*(1.+n_night),
                        c='C'+repr(n_night),
                        s=2, zorder=2,
                        label=night,
                        )

            ax2.errorbar(transmission_average['wave'],
                        transmission_average[night]['average_out']-spec_offset*(1.+n_night),
                        yerr=transmission_average[night]['average_out_err'],
                        color='C'+repr(n_night),
                        fmt='o', ms=1, zorder=1, alpha=0.25)

            ax2.scatter(transmission_average['wave'],
                        transmission_average[night]['average_out']-spec_offset*(1.+n_night),
                        c='C'+repr(n_night),
                        s=2, zorder=2)

        #ax1.set_ylim(0.95-spec_offset*(1.+n_night), 1.05)
        #ax1.set_xlim(config_in['master-out']['wavelength_range'][0], config_in['master-out']['wavelength_range'][1])
        try:
            ax1.set_xlim(lines_dict['plot_range'][0], lines_dict['plot_range'][1])
        except:
            ax1.set_xlim(lines_dict['range'][0], lines_dict['range'][1])

        ax1.set_ylim(0.985, 1.01)
        ax2.set_xlabel('$\lambda$ [$\AA$]')
        ax1.legend(loc=3)
        ax1.set_title('Average in-transit transmission spectrum in {0:s}'.format(reference))
        ax2.set_title('Average out-transit transmission spectrum in {0:s}'.format(reference))

        output_file = get_filename(filename_rad + '_binned', config_in['output'], night='', lines=lines_label, extension='.pdf', it_string=it_string)
        plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
        if interactive_plots:
            plt.show()
        plt.close()
