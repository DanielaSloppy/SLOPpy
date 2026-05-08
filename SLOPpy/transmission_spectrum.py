from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.math_functions import *
from SLOPpy.transmission_spectrum_preparation import compute_transmission_spectrum_preparation
from scipy.signal import savgol_filter

__all__ = ['compute_transmission_spectrum',
           'plot_transmission_spectrum']


subroutine_name = 'transmission_spectrum'
sampler_name = 'emcee'


def compute_transmission_spectrum(config_in, lines_label, reference='planetRF', night_input='', preparation_only=False, pca_iteration=-1):

    results_list_default = ['user',
                            'mcmc_night_MED',
                            'mcmc_night_MAP',
                            'mcmc_global_MED',
                            'mcmc_global_MAP']
    # compute_transmission_spectrum_preparation(config_in)

    night_dict = from_config_get_nights(config_in)
    ### transmission_dict = from_config_get_transmission(config_in)
    pca_parameters = from_config_get_pca_parameters(config_in)

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    norm_dict = lines_dict.get('normalization', {})
    norm_pams = {}
    norm_pams['normalize_transmission'] = norm_dict.get('normalize_transmission', True)
    norm_pams['normalization_model'] = norm_dict.get('normalization_model', 'polynomial')

    """ Normalization parameters for polynomial model"""
    norm_pams['model_poly_degree'] = norm_dict.get('model_poly_degree', 2)
    norm_pams['spectra_poly_degree'] = norm_dict.get('spectra_poly_degree', 2)
    norm_pams['lower_threshold'] = norm_dict.get('lower_threshold', 0.950)
    norm_pams['percentile_selection'] = norm_dict.get('percentile_selection', 10)

    """ Normalization parameters using Savitzky-Golay filter"""
    norm_pams['window_length'] = norm_dict.get('window_length', 101)
    norm_pams['polyorder'] = norm_dict.get('polyorder', 3)
    norm_pams['mode'] = norm_dict.get('mode', 'nearest')
    norm_pams['cval'] = norm_dict.get('cval', 1.0)

    shared_data = load_from_cpickle('shared', config_in['output'])

    clv_rm_correction = lines_dict.get('clv_rm_correction', True)

    """ Using the line-specific range to define the transmission spectrum region """

    shared_selection = (shared_data['coadd']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['coadd']['wave'] < lines_dict['range'][1])
    binned_selection = (shared_data['binned']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['binned']['wave'] < lines_dict['range'][1])

    transmission_template = {
        'subroutine': subroutine_name,
        'range': lines_dict['range'],
        'wave': shared_data['coadd']['wave'][shared_selection],
        'step': shared_data['coadd']['step'][shared_selection],
        'size': int(np.sum(shared_selection)),
        'binned_wave': shared_data['binned']['wave'][binned_selection],
        'binned_step': shared_data['binned']['step'][binned_selection],
        'binned_size': int(np.sum(binned_selection))
    }

    for night in night_dict:

        print()
        print("Running {0:45s} for  {1:20s}   Night:{2:15s}  ".format(subroutine_name, lines_label, night))

        preparation_input = load_from_cpickle('transmission_preparation', config_in['output'], night)

        if preparation_input.get('pca_output', False):
            if pca_iteration >= 0:
                it_string = str(pca_iteration).zfill(2)
            else:
                it_string = str(pca_parameters.get('ref_iteration', 3)).zfill(2)
            preparation = preparation_input[it_string]
        else:
            preparation = preparation_input
            it_string = ''

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        results_list = results_list_default.copy()

        binned_mcmc_night = check_existence_cpickle(
            'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label, it_string)
        binned_mcmc_global = check_existence_cpickle(
            'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], lines=lines_label, it_string=it_string)

        mcmc_night = check_existence_cpickle(
            'transmission_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label, it_string=it_string)
        mcmc_global = check_existence_cpickle(
            'transmission_mcmc_'+sampler_name+'_results', config_in['output'], lines=lines_label, it_string=it_string)

        if mcmc_night and mcmc_global:
            mcmc_results_night = load_from_cpickle(
                'transmission_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label, it_string=it_string)

            mcmc_results_global = load_from_cpickle(
                'transmission_mcmc_'+sampler_name+'_results', config_in['output'], lines=lines_label, it_string=it_string)

            print('  Observational parameters from MCMC fit of unbinned data and configuration file')

        elif binned_mcmc_night and binned_mcmc_global:
            mcmc_results_night = load_from_cpickle(
                'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label, it_string=it_string)

            mcmc_results_global = load_from_cpickle(
                'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], lines=lines_label, it_string=it_string)

            print('  Observational parameters from MCMC fit of binned data and configuration file')

        else:
            print('  Observational parameters from configuration file')
            results_list = ['user']

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        if clv_rm_correction:
            try:
                clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night, lines_label)
            except (FileNotFoundError, IOError):
                clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night)

        for results_selection in results_list:

            try:
                transmission = load_from_cpickle(subroutine_name+'_'+reference + '_' +
                                                 results_selection, config_in['output'], night, lines_label, it_string=it_string)
                print("{0:45s} Night:{1:15s}   {2:s}   {3:s}   {4:s}".format(
                    subroutine_name, night, lines_label, results_selection, 'Retrieved'))
                continue
            except (FileNotFoundError, IOError):
                print("{0:45s} Night:{1:15s}   {2:s}   {3:s}   {4:s}".format(
                    subroutine_name, night, lines_label, results_selection, 'Computing'))

            transmission = transmission_template.copy()

            if len(it_string) > 0:
                transmission['pca_output'] = True
            else:
                transmission['pca_output'] = False

            print_warning = True

            for obs in lists['observations']:

                """ we start from the e2ds file, after correction for blaze and
                    division by the master-out

                Observation data:
                    wave: input_data[obs]['wave']
                    step: input_data[obs]['step']
                    flux: preparation[obs]['deblazed']
                    ferr: preparation[obs]['deblazed_err']
                """

                transmission[obs] = {}
                transmission[obs] = {
                    'BJD': input_data[obs]['BJD'],
                    'AIRMASS': input_data[obs]['AIRMASS']
                }

                """ Shift into planetary reference system is the default
                choice"""
                if results_selection == 'user':
                    planet_R_factor = observational_pams.get('Rp_factor', 1.00000)
                    if reference in ['observer', 'observerRF', 'ORF']:
                        rv_shift = 0.000
                        rv_shift_clv = -observational_pams[obs]['rv_shift_ORF2SRF']
                    elif reference in ['stellar', 'stellarRF', 'SRF']:
                        rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                        rv_shift_clv = 0.0000
                    else:
                        rv_shift = observational_pams[obs]['rv_shift_ORF2PRF']
                        rv_shift_clv = observational_pams[obs]['rv_shift_SRF2PRF']

                elif results_selection == 'mcmc_night_MED':
                    planet_R_factor = mcmc_results_night['results']['planet_R']
                    if reference in ['observer', 'observerRF', 'ORF']:
                        rv_shift = 0.000
                        #rv_shift_clv = mcmc_results_night['results']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift_clv = -observational_pams[obs]['rv_shift_ORF2SRF']
                    elif reference in ['stellar', 'stellarRF', 'SRF']:
                        #rv_shift = mcmc_results_night['results']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                        rv_shift_clv = 0.0000
                    else:
                        rv_shift = mcmc_results_night['results']['observational_pams'][obs]['rv_shift_ORF2PRF']
                        rv_shift_clv = mcmc_results_night['results']['observational_pams'][obs]['rv_shift_SRF2PRF']

                elif results_selection == 'mcmc_night_MAP':
                    planet_R_factor = mcmc_results_night['results_MAP']['planet_R']
                    if reference in ['observer', 'observerRF', 'ORF']:
                        rv_shift = 0.000
                        #rv_shift_clv = mcmc_results_night['results_MAP']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = -observational_pams[obs]['rv_shift_ORF2SRF']
                    elif reference in ['stellar', 'stellarRF', 'SRF']:
                        #rv_shift = mcmc_results_night['results_MAP']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                        rv_shift_clv = 0.0000
                    else:
                        rv_shift = mcmc_results_night['results_MAP']['observational_pams'][obs]['rv_shift_ORF2PRF']
                        rv_shift_clv = mcmc_results_night['results_MAP']['observational_pams'][obs]['rv_shift_SRF2PRF']

                elif results_selection == 'mcmc_global_MED':
                    planet_R_factor = mcmc_results_global['results']['planet_R']
                    if reference in ['observer', 'observerRF', 'ORF']:
                        rv_shift = 0.000
                        #rv_shift_clv = mcmc_results_global['results']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = -observational_pams[obs]['rv_shift_ORF2SRF']
                    elif reference in ['stellar', 'stellarRF', 'SRF']:
                        #rv_shift = mcmc_results_global['results']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                        rv_shift_clv = 0.0000
                    else:
                        rv_shift = mcmc_results_global['results']['observational_pams'][obs]['rv_shift_ORF2PRF']
                        rv_shift_clv = mcmc_results_global['results']['observational_pams'][obs]['rv_shift_SRF2PRF']

                elif results_selection == 'mcmc_global_MAP':
                    planet_R_factor = mcmc_results_global['results_MAP']['planet_R']
                    if reference in ['observer', 'observerRF', 'ORF']:
                        rv_shift = 0.000
                        #rv_shift_clv = mcmc_results_global['results_MAP']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = -observational_pams[obs]['rv_shift_ORF2SRF']
                    elif reference in ['stellar', 'stellarRF', 'SRF']:
                        #rv_shift = mcmc_results_global['results_MAP']['observational_pams'][obs]['rv_shift_ORF2SRF']
                        rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                        rv_shift_clv = 0.0000
                    else:
                        rv_shift = mcmc_results_global['results_MAP']['observational_pams'][obs]['rv_shift_ORF2PRF']
                        rv_shift_clv = mcmc_results_global['results_MAP']['observational_pams'][obs]['rv_shift_SRF2PRF']

                """ Step 2): rebin the 2D ratio spectra to 1D """

                if transmission['pca_output']:

                    transmission[obs]['rebinned'] = \
                        rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio'],
                                    np.ones_like(calib_data['blaze']),
                                    transmission['wave'],
                                    transmission['step'],
                                    preserve_flux=False,
                                    rv_shift=rv_shift)

                    transmission[obs]['rebinned_err'] = \
                        rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio_err'],
                                    np.ones_like(calib_data['blaze']),
                                    transmission['wave'],
                                    transmission['step'],
                                    rv_shift=rv_shift,
                                    preserve_flux=False,
                                    is_error=True)
                else:
                    # preserve_flux is intentionally forced False for the ratio:
                    # ratio is dimensionless, and preserve_flux=True embeds the
                    # variable ESPRESSO pixel width into the rebinned ratio,
                    # producing a spurious ~10% slope that corrupts normalization.
                    transmission[obs]['rebinned'] = \
                        rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio'],
                                    calib_data['blaze'],
                                    transmission['wave'],
                                    transmission['step'],
                                    preserve_flux=False,
                                    rv_shift=rv_shift)

                    transmission[obs]['rebinned_err'] = \
                        rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio_err'],
                                    calib_data['blaze'],
                                    transmission['wave'],
                                    transmission['step'],
                                    preserve_flux=False,
                                    rv_shift=rv_shift,
                                    is_error=True)

                    ### Small border bugfix
                    if transmission[obs]['rebinned_err'][0] ==0:
                        transmission[obs]['rebinned'][0] = transmission[obs]['rebinned'][1]
                        transmission[obs]['rebinned_err'][0] = transmission[obs]['rebinned_err'][1]

                    if transmission[obs]['rebinned_err'][-1] ==0:

                        transmission[obs]['rebinned'][-1] = transmission[obs]['rebinned'][-2]
                        transmission[obs]['rebinned_err'][-1] = transmission[obs]['rebinned_err'][-2]


                    #import matplotlib.pyplot as plt
                    #plt.scatter(transmission['wave'], transmission[obs]['corrected'])
                    #plt.plot(transmission['wave'], transmission[obs]['continuum'])
                    #plt.scatter(transmission['wave'][selection], transmission[obs]['corrected'][selection], c='r')
                    #plt.plot(input_data[obs]['wave'][0,:], preparation[obs]['ratio_err'][0,:])
                    #plt.scatter(transmission['wave'], transmission[obs]['rebinned_err'], c='b')
                    #plt.axhline(0.0000, c='C2')

                    #plt.show()
                    #quit()

                #import matplotlib.pyplot as plt
                #plt.scatter(input_data[obs]['wave'], preparation[obs]['ratio'], s=2)
                #plt.xlim(lines_dict['range'][0], lines_dict['range'][1])
                # plt.show()

                if clv_rm_correction:

                    """" CLV + RM computation in the planetary reference frame """
                    transmission[obs]['clv_model_stellarRF'] = interpolate1d_grid_nocheck(planet_R_factor,
                                                                                          clv_rm_models['common']['radius_grid'],
                                                                                          clv_rm_models[obs]['clv_rm_model_convolved_normalized'])

                    transmission[obs]['clv_model_rebinned'] = \
                        rebin_1d_to_1d(clv_rm_models['common']['wave'],
                                       clv_rm_models['common']['step'],
                                       transmission[obs]['clv_model_stellarRF'],
                                       transmission['wave'],
                                       transmission['step'],
                                       preserve_flux=False,
                                       rv_shift=rv_shift_clv)

                    #import matplotlib.pyplot as plt
                    #print(obs, planet_R_factor)
                    #plt.plot(clv_rm_models['common']['wave'], transmission[obs]['clv_model_stellarRF'], zorder=100, c='C2')
                    #plt.scatter(transmission['wave'], transmission[obs]['clv_model_rebinned'], s=2)
                    # plt.show()
                    transmission[obs]['corrected'] = transmission[obs]['rebinned'] / \
                        transmission[obs]['clv_model_rebinned']
                    transmission[obs]['corrected_err'] = transmission[obs]['rebinned_err'] / \
                        transmission[obs]['clv_model_rebinned']

                else:
                    transmission[obs]['clv_model_rebinned'] = np.ones(transmission['size'])
                    transmission[obs]['corrected'] = transmission[obs]['rebinned']
                    transmission[obs]['corrected_err'] = transmission[obs]['rebinned_err']
                    if print_warning:
                        print('   *** No CLV correction')


                if norm_pams['normalize_transmission'] and norm_pams['normalization_model'] == 'polynomial':

                    """ Continuum normalization preparatory steps:
                        1) exclusion of regions with lines of interes
                        2) exclusion of regions with stellar lines
                        3) Polynomial fit of selected regions
                        Boolean array initialized to all True values
                    """
                    transmission[obs]['line_exclusion'] = (transmission['wave'] > 0.)

                    """ Continuum normalization:
                        1) exclusion of regions with transmission lines under study, now
                        in the RF of the lines
                    """
                    for line_key, line_val in lines_dict['lines'].items():
                        transmission[obs]['line_exclusion'] = transmission[obs]['line_exclusion'] & (
                            np.abs(transmission['wave']-line_val) > 3.)

                    """ Continuum normalization:
                        2) exclusion of regions with planetary lines, taking into account the planetary RV semi-amplitude
                    """

                    if clv_rm_correction:
                        stellar_spectrum_rebinned = rebin_1d_to_1d(clv_rm_models['common']['wave'],
                                                                clv_rm_models['common']['step'],
                                                                clv_rm_models['common']['norm_convolved'],
                                                                transmission['wave'],
                                                                transmission['step'],
                                                                rv_shift=rv_shift_clv,
                                                                preserve_flux=False)

                        stellar_spectrum_derivative = first_derivative(transmission['wave'], stellar_spectrum_rebinned)

                        missing_model = (np.abs(stellar_spectrum_rebinned) < 0.0001)

                        cont_10perc = np.percentile(np.abs(stellar_spectrum_derivative), norm_pams['percentile_selection'])

                        #transmission[obs]['line_exclusion'] = transmission[obs]['line_exclusion'] \
                        #    & (np.abs(stellar_spectrum_derivative) < cont_10perc) \
                        #    & (stellar_spectrum_rebinned > norm_pams['lower_threshold'])

                        line_exclusion = transmission[obs]['line_exclusion'] \
                            & (np.abs(stellar_spectrum_derivative) < cont_10perc) \
                            & (stellar_spectrum_rebinned > norm_pams['lower_threshold'])

                        if np.sum(line_exclusion) < len(line_exclusion)/200:
                            transmission[obs]['line_exclusion'] = transmission[obs]['line_exclusion'] \
                                & ( missing_model | ((np.abs(stellar_spectrum_derivative) < cont_10perc) \
                                & (stellar_spectrum_rebinned > norm_pams['lower_threshold'])))
                        else:
                            transmission[obs]['line_exclusion'] = line_exclusion



                    elif print_warning:
                        print("   No stellar synthetic spectrum from CLV models")
                        print("   some stellar lines may be included in transmission normalization  ")
                        print_warning = False

                    """ Continuum normalization:
                        3) Polynomial fit, everything is hard coded now but personalized
                        options can be implemented easily in the yaml file
                    """

                    selection = transmission[obs]['line_exclusion'] & (
                        transmission[obs]['corrected'] > np.std(transmission[obs]['corrected']))
                    transmission[obs]['continuum_coeff'] = \
                        np.polynomial.chebyshev.chebfit(transmission['wave'][selection],
                                                        transmission[obs]['corrected'][selection],
                                                        norm_pams['spectra_poly_degree'])
                    transmission[obs]['continuum'] = np.polynomial.chebyshev.chebval(
                        transmission['wave'], transmission[obs]['continuum_coeff'])

                    transmission[obs]['normalized'] = transmission[obs]['corrected'] / transmission[obs]['continuum']
                    transmission[obs]['normalized_err'] = transmission[obs]['corrected_err'] / \
                        transmission[obs]['continuum']

                    #import matplotlib.pyplot as plt
                    #plt.scatter(transmission['wave'], transmission[obs]['corrected'])
                    #plt.plot(transmission['wave'], transmission[obs]['continuum'])
                    #plt.scatter(transmission['wave'][selection], transmission[obs]['corrected'][selection], c='r')
                    #plt.scatter(transmission['wave'], transmission[obs]['corrected_err']+0.05, c='b')
                    #plt.scatter(transmission['wave'], transmission[obs]['normalized_err'], c='r')
                    #plt.show()
                    #quit()

                    transmission[obs]['continuum_uncorrected_coeff'] = \
                        np.polynomial.chebyshev.chebfit(transmission['wave'][selection],
                                                        transmission[obs]['rebinned'][selection],
                                                        norm_pams['spectra_poly_degree'])
                    transmission[obs]['continuum_uncorrected'] = np.polynomial.chebyshev.chebval(
                        transmission['wave'], transmission[obs]['continuum_uncorrected_coeff'])
                    transmission[obs]['normalized_uncorrected'] = transmission[obs]['rebinned'] / \
                        transmission[obs]['continuum_uncorrected']
                    transmission[obs]['normalized_uncorrected_err'] = transmission[obs]['rebinned_err'] / \
                        transmission[obs]['continuum_uncorrected']

                elif norm_pams['normalize_transmission']  and (
                    norm_pams['normalization_model'] == 'savgol'
                    or norm_pams['normalization_model'] == 'savitzky-golay'):

                    print(' ', obs, ' normalization using Savitzky-Golay filter')

                    transmission[obs]['continuum_coeff'] = None
                    transmission[obs]['continuum_uncorrected_coeff'] = None


                    transmission[obs]['continuum'] = savgol_filter(transmission[obs]['corrected'],
                                                            window_length=norm_pams['window_length'],
                                                            polyorder=norm_pams['polyorder'],
                                                            mode=norm_pams['mode'],
                                                            cval=norm_pams['cval'])

                    transmission[obs]['normalized'] = transmission[obs]['corrected'] / transmission[obs]['continuum']
                    transmission[obs]['normalized_err'] = transmission[obs]['corrected_err'] / \
                        transmission[obs]['continuum']

                    transmission[obs]['continuum_uncorrected'] = savgol_filter(transmission[obs]['rebinned'],
                                                            window_length=norm_pams['window_length'],
                                                            polyorder=norm_pams['polyorder'],
                                                            mode=norm_pams['mode'],
                                                            cval=norm_pams['cval'])

                    transmission[obs]['normalized_uncorrected'] = transmission[obs]['rebinned'] / transmission[obs]['continuum_uncorrected']
                    transmission[obs]['normalized_uncorrected_err'] = transmission[obs]['rebinned_err'] / \
                        transmission[obs]['continuum_uncorrected']

                else:
                    transmission[obs]['continuum_coeff'] = None
                    transmission[obs]['continuum'] = np.ones_like(transmission['wave'])
                    transmission[obs]['normalized'] = transmission[obs]['corrected'].copy()
                    transmission[obs]['normalized_err'] = transmission[obs]['corrected_err'].copy()

                    #import matplotlib.pyplot as plt
                    #plt.scatter(transmission['wave'], transmission[obs]['corrected'])
                    #plt.plot(transmission['wave'], transmission[obs]['continuum'])
                    #plt.scatter(transmission['wave'][selection], transmission[obs]['corrected'][selection], c='r')
                    # plt.show()

                    transmission[obs]['continuum_uncorrected_coeff'] = None
                    transmission[obs]['continuum_uncorrected'] = np.ones_like(transmission['wave'])
                    transmission[obs]['normalized_uncorrected'] = transmission[obs]['rebinned'].copy()
                    transmission[obs]['normalized_uncorrected_err'] = transmission[obs]['rebinned_err'].copy()
                    print_warning = False

            transm_average = np.zeros([len(lists['transit_full']), transmission['size']])
            weights_average = np.zeros([len(lists['transit_full']), transmission['size']])
            clvrm_average = np.zeros([len(lists['transit_full']), transmission['size']])
            uncorr_average = np.zeros([len(lists['transit_full']), transmission['size']])

            for i, obs in enumerate(lists['transit_full']):
                transm_average[i, :] = transmission[obs]['normalized'][:]
                weights_average[i, :] = 1./(transmission[obs]['normalized_err']**2.)
                clvrm_average[i, :] = transmission[obs]['clv_model_rebinned'][:]
                uncorr_average[i, :] = transmission[obs]['normalized_uncorrected'][:]

            transmission['average'], transmission['sum_weights'] = np.average(
                transm_average, axis=0, weights=weights_average, returned=True)
            transmission['average_err'] = 1. / np.sqrt(transmission['sum_weights'])

            transmission['average_clv_model'], _ = np.average(
                clvrm_average, axis=0, weights=weights_average, returned=True)

            transmission['average_uncorrected'], _ = np.average(
                uncorr_average, axis=0, weights=weights_average, returned=True)

            transmission['binned'] = \
                rebin_1d_to_1d(transmission['wave'],
                               transmission['step'],
                               transmission['average'],
                               transmission['binned_wave'],
                               transmission['binned_step'],
                               preserve_flux=False)
            transmission['binned_err'] = \
                rebin_1d_to_1d(transmission['wave'],
                               transmission['step'],
                               transmission['average_err'],
                               transmission['binned_wave'],
                               transmission['binned_step'],
                               preserve_flux=False,
                               is_error=True)

            transmission['binned_clv_model'] = \
                rebin_1d_to_1d(transmission['wave'],
                               transmission['step'],
                               transmission['average_clv_model'],
                               transmission['binned_wave'],
                               transmission['binned_step'],
                               preserve_flux=False)

            transmission['binned_uncorrected'] = \
                rebin_1d_to_1d(transmission['wave'],
                               transmission['step'],
                               transmission['average_uncorrected'],
                               transmission['binned_wave'],
                               transmission['binned_step'],
                               preserve_flux=False)

            transm_average = np.zeros([len(lists['transit_out']), transmission['size']])
            weights_average = np.zeros([len(lists['transit_out']), transmission['size']])

            for i, obs in enumerate(lists['transit_out']):
                transm_average[i, :] = transmission[obs]['normalized'][:]
                weights_average[i, :] = 1./(transmission[obs]['normalized_err']**2.)

            transmission['average_out'], transmission['sum_weights_out'] = np.average(
                transm_average, axis=0, weights=weights_average, returned=True)
            transmission['average_out_err'] = 1./np.sqrt(transmission['sum_weights_out'])

            transmission['binned_out'] = \
                rebin_1d_to_1d(transmission['wave'],
                               transmission['step'],
                               transmission['average_out'],
                               transmission['binned_wave'],
                               transmission['binned_step'],
                               preserve_flux=False)
            transmission['binned_out_err'] = \
                rebin_1d_to_1d(transmission['wave'],
                               transmission['step'],
                               transmission['average_out_err'],
                               transmission['binned_wave'],
                               transmission['binned_step'],
                               preserve_flux=False,
                               is_error=True)

            #save_to_cpickle('transmission_'+reference+'_processed', processed, config_in['output'], night)
            save_to_cpickle(subroutine_name + '_' + reference + '_' + results_selection,
                            transmission, config_in['output'], night, lines_label, it_string)

            # Forcing memory deallocation
            transmission = None

        # Forcing memory deallocation
        clv_rm_models = None


def plot_transmission_spectrum(config_in, lines_label, night_input='', results_input='', reference='planetRF', pca_iteration=-1):

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    if results_input == '':
        results_list = ['user',
                        'mcmc_night_MED',
                        'mcmc_night_MAP',
                        'mcmc_global_MED',
                        'mcmc_global_MAP']
    else:
        results_list = np.atleast_1d(results_input)

    clv_rm_correction = lines_dict.get('clv_rm_correction', True)

    os.system('mkdir -p plots')

    interactive_plots = from_config_get_interactive_plots(config_in)

    for night in night_list:

        # Workaround to check if the transmission spectrum has been obtained through PCA iterations
        preparation_input = load_from_cpickle('transmission_preparation', config_in['output'], night)

        if preparation_input.get('pca_output', False):
            if pca_iteration >= 0:
                it_string = str(pca_iteration).zfill(2)
            else:
                it_string = str(preparation_input.get('ref_iteration', 0)).zfill(2)
        else:
            it_string = ''
        preparation_input = None

        if clv_rm_correction:
            try:
                clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night, lines_label)
            except (FileNotFoundError, IOError):
                clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night)

        for results_selection in results_list:

            filename_rad = subroutine_name + '_'+reference+'_'+results_selection

            """ Retrieving the list of observations"""
            lists = load_from_cpickle('lists', config_in['output'], night)

            """ Retrieving the analysis"""
            try:
                #processed = load_from_cpickle('transmission_'+reference+'_processed', config_in['output'], night)
                transmission = load_from_cpickle(filename_rad, config_in['output'], night, lines_label, it_string)
            except (FileNotFoundError, IOError):
                print()
                print("No transmission spectrum in {0:s}, no plots".format(reference))
                continue

            """ Creation of the color array, based on the BJD of the observations
            """
            bjd = []
            am = []

            for obs in lists['observations']:
                bjd.append(transmission[obs]['BJD'] - 2450000.0)
                am.append(transmission[obs]['AIRMASS'])

            color_cmap = plt.cm.viridis
            color_norm = plt.Normalize(vmin=bjd[0], vmax=bjd[-1])
            colors = color_cmap(color_norm(np.asarray(bjd)))

            fig = plt.figure(figsize=(12, 6))

            gs = GridSpec(2, 2, width_ratios=[50, 1])
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
            cbax1 = plt.subplot(gs[:, 1])

            # commented out because the plot was too cumbersome
            for obs in lists['transit_full']:

                color = [color_cmap(color_norm(transmission[obs]['BJD'] - 2450000.0))[:-1]]

                ax1.scatter(transmission['wave'],
                            transmission[obs]['normalized'],
                            c=color, s=1, zorder=3, alpha=0.25)

            for obs in lists['transit_out']:

                color = [color_cmap(color_norm(transmission[obs]['BJD'] - 2450000.0))[:-1]]

                ax2.scatter(transmission['wave'],
                            transmission[obs]['normalized'],
                            c=color, s=1, zorder=3, alpha=0.25)

            ax1.set_ylim(0.925, 1.075)
            ax2.set_xlabel('$\lambda$ [$\AA$]')
            ax2.legend(loc=3)
            ax1.set_title('Lines: {0:s} Night: {1:s} \n In-transit transmission spectrum in {2:s} \n Solution {3:s}'.format(
                lines_label, night, reference, results_selection))
            ax2.set_title('Out-transit transmission spectrum in {0:s}'.format(reference))
            try:
                ax1.set_xlim(lines_dict['plot_range'][0], lines_dict['plot_range'][1])
            except:
                ax1.set_xlim(lines_dict['range'][0], lines_dict['range'][1])

            sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
            sm.set_array([])  # You have to set a dummy-array for this to work...
            cbar = plt.colorbar(sm, cax=cbax1)
            cbar.set_label('BJD - 2450000.0')
            fig.subplots_adjust(wspace=0.05, hspace=0.4)

            output_file = get_filename(filename_rad + '_observations',
                                       config_in['output'], night, lines_label, it_string, extension='.pdf')
            plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
            if interactive_plots:
                plt.show()
            plt.close()

            fig = plt.figure(figsize=(12, 6))

            gs = GridSpec(2, 2, width_ratios=[50, 1])
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
            cbax1 = plt.subplot(gs[:, 1])

            try:
                master_out = load_from_cpickle('master_out', config_in['output'], night)
                ax2.plot(master_out['wave'],
                         master_out['rescaled']-0.06,
                         color='k', zorder=10, label='master-out')
            except (FileNotFoundError, IOError):
                pass

            try:
                telluric = load_from_cpickle('telluric', config_in['output'], night)
                ax2.plot(telluric['template']['input']['wave'],
                         telluric['template']['input']['flux'] - 0.06,
                         color='C1', zorder=10, label='telluric')
                ax2.plot(telluric['template']['input']['wave'],
                         (telluric['template']['input']['flux']-1.)*10. + 1. - 0.06,
                         color='C2', alpha=0.5, zorder=9, label='telluric (x10)')
            except (FileNotFoundError, IOError, KeyError):
                pass

            #master_out = load_from_cpickle('master_out', config_in['output'], night)
            # ax1.errorbar(master_out['wave'],
            #            master_out['rescaled'],
            #            yerr=master_out['rescaled_err'],
            #            fmt='.', c='C0', label='master-out ' + night)

            ax1.errorbar(transmission['wave'],
                         transmission['average'],
                         yerr=transmission['average_err'],
                         fmt='ko', ms=1, zorder=5, alpha=0.25)

            ax1.errorbar(transmission['binned_wave'],
                         transmission['binned'],
                         yerr=transmission['binned_err'],
                         fmt='ro', ms=4, lw=2,  zorder=10)

            ax2.errorbar(transmission['wave'],
                         transmission['average_out'],
                         yerr=transmission['average_out_err'],
                         fmt='ko', ms=1, zorder=5, alpha=0.25, label='average')

            ax2.errorbar(transmission['binned_wave'],
                         transmission['binned_out'],
                         yerr=transmission['binned_out_err'],
                         fmt='ro', ms=4, lw=2,  zorder=10, label='binned average')

            ax1.set_ylim(0.99, 1.01)
            ax2.set_xlabel('$\lambda$ [$\AA$]')
            ax2.legend(loc=3)
            ax1.set_title('Lines: {0:s} Night: {1:s} \n In-transit transmission spectrum in {2:s} \n Solution {3:s}'.format(
                lines_label, night, reference, results_selection))
            ax2.set_title('Out-transit transmission spectrum in {0:s}'.format(reference))

            sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
            sm.set_array([])  # You have to set a dummy-array for this to work...
            cbar = plt.colorbar(sm, cax=cbax1)
            cbar.set_label('BJD - 2450000.0')
            fig.subplots_adjust(wspace=0.05, hspace=0.4)

            try:
                ax1.set_xlim(lines_dict['plot_range'][0], lines_dict['plot_range'][1])
            except:
                ax1.set_xlim(lines_dict['range'][0], lines_dict['range'][1])

            #ax1.set_xlim(config_in['master-out']['wavelength_range'][0], config_in['master-out']['wavelength_range'][1])
            output_file = get_filename(filename_rad + '_binned',
                                       config_in['output'], night, lines_label, it_string, extension='.pdf')
            plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
            if interactive_plots:
                plt.show()
            plt.close()

            if not clv_rm_correction:
                continue

            fig = plt.figure(figsize=(12, 6))

            gs = GridSpec(2, 2, width_ratios=[50, 1])
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
            cbax1 = plt.subplot(gs[:, 1])

            # commented out because the plot was too cumbersome
            for obs in lists['transit_full']:

                color = [color_cmap(color_norm(transmission[obs]['BJD'] - 2450000.0))[:-1]]

                ax1.plot(clv_rm_models['common']['wave'],
                         transmission[obs]['clv_model_stellarRF'],
                         zorder=3, alpha=0.25)

                ax1.scatter(transmission['wave'],
                            transmission[obs]['clv_model_rebinned'],
                            c=color, s=1, zorder=10, alpha=0.5)

            for obs in lists['transit_out']:

                color = [color_cmap(color_norm(transmission[obs]['BJD'] - 2450000.0))[:-1]]

                ax2.plot(clv_rm_models['common']['wave'],
                         transmission[obs]['clv_model_stellarRF'],
                         zorder=3, alpha=0.25)

                ax2.scatter(transmission['wave'],
                            transmission[obs]['clv_model_rebinned'],
                            c=color, s=1, zorder=10, alpha=0.5)

            ax2.set_xlabel('$\lambda$ [$\AA$]')
            ax2.legend(loc=3)
            ax1.set_title('Lines: {0:s} Night: {1:s} \n CLV-RM correction in {2:s} \n Solution {3:s}'.format(
                lines_label, night, reference, results_selection))
            ax2.set_title('Out-transit transmission spectrum in {0:s}'.format(reference))
            try:
                ax1.set_xlim(lines_dict['plot_range'][0], lines_dict['plot_range'][1])
            except:
                ax1.set_xlim(lines_dict['range'][0], lines_dict['range'][1])

            sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
            sm.set_array([])  # You have to set a dummy-array for this to work...
            cbar = plt.colorbar(sm, cax=cbax1)
            cbar.set_label('BJD - 2450000.0')
            fig.subplots_adjust(wspace=0.05, hspace=0.4)

            output_file = get_filename(filename_rad + '_clv_rm_models',
                                       config_in['output'], night, lines_label, it_string, extension='.pdf')
            plt.savefig('plots/'+output_file, bbox_inches='tight', dpi=300)
            if interactive_plots:
                plt.show()
            plt.close()
