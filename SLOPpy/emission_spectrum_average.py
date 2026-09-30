from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.emission_spectrum import emission_template_from_range, average_and_bin, emission_ylim

__all__ = ['compute_emission_spectrum_average',
           'plot_emission_spectrum_average',
           'compute_emission_spectrum_average_planetRF',
           'plot_emission_spectrum_average_planetRF',
           'compute_emission_spectrum_average_stellarRF',
           'plot_emission_spectrum_average_stellarRF',
           'compute_emission_spectrum_average_observerRF',
           'plot_emission_spectrum_average_observerRF']


subroutine_name = 'emission_spectrum_average'
pick_files = 'emission_spectrum'


def compute_emission_spectrum_average_planetRF(config_in, lines_label):
    compute_emission_spectrum_average(config_in, lines_label, reference='planetRF')


def plot_emission_spectrum_average_planetRF(config_in, lines_label):
    plot_emission_spectrum_average(config_in, lines_label, reference='planetRF')


def compute_emission_spectrum_average_stellarRF(config_in, lines_label):
    compute_emission_spectrum_average(config_in, lines_label, reference='stellarRF')


def plot_emission_spectrum_average_stellarRF(config_in, lines_label):
    plot_emission_spectrum_average(config_in, lines_label, reference='stellarRF')


def compute_emission_spectrum_average_observerRF(config_in, lines_label):
    compute_emission_spectrum_average(config_in, lines_label, reference='observerRF')


def plot_emission_spectrum_average_observerRF(config_in, lines_label):
    plot_emission_spectrum_average(config_in, lines_label, reference='observerRF')


def compute_emission_spectrum_average(config_in, lines_label, reference='planetRF'):
    """
    Average emission spectrum of all the eclipse nights. The spectra of the individual
    observations of every night are combined with a weighted average.
    The average is always recomputed, as it only requires the emission spectra of the nights
    """

    night_dict = from_config_get_nights(config_in, phase='eclipse')

    if len(night_dict) == 0:
        print("{0:45s}    {1:s}".format(subroutine_name + '_' + reference, 'No eclipse nights, skipped'))
        return

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    shared_data = load_from_cpickle('shared', config_in['output'])
    emission_average = emission_template_from_range(shared_data, subroutine_name + '_' + reference,
                                                    lines_dict['range'])
    emission_average['reference'] = reference
    emission_average['nights'] = list(night_dict)

    spectra = {'': [], '_out': []}
    spectra_err = {'': [], '_out': []}

    for night in night_dict:

        lists = load_from_cpickle('lists', config_in['output'], night)

        try:
            emission_average[night] = load_from_cpickle(pick_files + '_' + reference, config_in['output'],
                                                        night, lines_label)
        except (FileNotFoundError, IOError):
            print("{0:45s} Night:{1:15s}   {2:s}   {3:s}".format(
                subroutine_name + '_' + reference, night, lines_label, 'Missing emission spectrum, skipped'))
            return

        for prefix, list_name in [['', 'transit_full'], ['_out', 'transit_out']]:
            for obs in lists[list_name]:
                spectra[prefix].append(emission_average[night][obs]['rebinned'])
                spectra_err[prefix].append(emission_average[night][obs]['rebinned_err'])

    print("{0:45s}    {1:s}   {2:s}".format(subroutine_name + '_' + reference, lines_label, 'Computing'))

    for prefix in ['', '_out']:
        average_and_bin(emission_average, prefix, spectra[prefix], spectra_err[prefix])

    save_to_cpickle(subroutine_name + '_' + reference, emission_average, config_in['output'], lines=lines_label)


def plot_emission_spectrum_average(config_in, lines_label, reference='planetRF'):

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]
    plot_range = lines_dict.get('plot_range', lines_dict['range'])

    os.system('mkdir -p plots')

    interactive_plots = from_config_get_interactive_plots(config_in)

    filename_rad = subroutine_name + '_' + reference

    try:
        emission_average = load_from_cpickle(filename_rad, config_in['output'], lines=lines_label)
        print("{0:45s}    {1:s}   {2:s}".format(filename_rad, lines_label, 'Plotting'))
    except (FileNotFoundError, IOError):
        print("{0:45s}                         {1:s}".format(filename_rad, 'Plot skipped'))
        return

    fig = plt.figure(figsize=(12, 9))

    gs = GridSpec(2, 1)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)

    ylim = emission_ylim(emission_average['average'], emission_average['average_err'])
    spec_offset = ylim[1]

    for ax, prefix in [[ax1, ''], [ax2, '_out']]:
        ax.errorbar(emission_average['wave'],
                    emission_average['average' + prefix],
                    yerr=emission_average['average' + prefix + '_err'],
                    fmt='ko', ms=1, zorder=10, alpha=0.10, label='average')

        ax.errorbar(emission_average['binned_wave'],
                    emission_average['binned' + prefix],
                    yerr=emission_average['binned' + prefix + '_err'],
                    fmt='ko', ms=3, zorder=20, label='binned')

        """ average of each night, shifted below the global average """
        for n_night, night in enumerate(emission_average['nights']):
            ax.errorbar(emission_average['wave'],
                        emission_average[night]['average' + prefix] - spec_offset*(1.+n_night),
                        yerr=emission_average[night]['average' + prefix + '_err'],
                        color='C'+repr(n_night),
                        fmt='o', ms=1, zorder=1, alpha=0.25)

            ax.scatter(emission_average['wave'],
                       emission_average[night]['average' + prefix] - spec_offset*(1.+n_night),
                       c='C'+repr(n_night),
                       s=2, zorder=2,
                       label=night)

    ax1.set_xlim(plot_range[0], plot_range[1])
    ax1.set_ylim(ylim[0] - spec_offset*len(emission_average['nights']), ylim[1])
    ax2.set_xlabel(r'$\lambda$ [$\AA$]')
    ax1.legend(loc=3)
    ax1.set_title('Average emission spectrum (full_transit list) in {0:s}'.format(reference))
    ax2.set_title('Average emission spectrum (out_transit list) in {0:s}'.format(reference))

    output_file = get_filename(filename_rad + '_binned', config_in['output'], night='', lines=lines_label,
                               extension='.pdf')
    plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
    if interactive_plots:
        plt.show()
    plt.close()
