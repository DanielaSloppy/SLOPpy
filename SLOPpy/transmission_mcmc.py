from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.constants import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.math_functions import *
from SLOPpy.subroutines.bayesian_emcee import *
# from SLOPpy.subroutines.rebin_subroutines import *
from scipy.signal import savgol_filter

__all__ = ['compute_transmission_mcmc', 'compute_transmission_mcmc_iterative']

subroutine_name = 'transmission_mcmc'


def compute_transmission_mcmc_iterative(config_in, lines_label):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0, pca_parameters.get('iterations', 5)):
        compute_transmission_mcmc(config_in, lines_label, reference='planetRF', pca_iteration=it)


def compute_transmission_mcmc(config_in, lines_label, reference='planetRF', pca_iteration=-1):

    night_dict = from_config_get_nights(config_in)
    planet_dict = from_config_get_planet(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    """ selection of those parameters that are specific of the spectral line(s)
        under analysis
    """

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

    sampler_pams = lines_dict['sampler_parameters']
    sampler_name = sampler_pams.get('sampler_name', 'emcee')

    # TODO reference as input parameter
    reference = 'planetRF'

    """
    - case  0: only one spectral line, default line parameters are contrast, FWHM, rv_shift
    - case  1: only one spectral line, no winds
    - case  2: only one spectral line, no planetary radius dependance
    - case  3: only one spectral line, no winds and no planetary radius dependance
    - case 10: more than one spectral lines, all line parameters are free and independent
    - case 11: more than one spectral lines, all lines are affected by the same wind
    - case 12: more than one spectral lines, all lines have same FWHM
    - case 13: more than one spectral lines, all lines are affected by the same wind and have same FWHM
    - case 14: more than one spectral lines, no winds
    - case 15: more than one spectral lines, no winds, all lines have same FWHM
    - case 20: more than one spectral lines, no Rp dependance, all line parameters are free and independent
    - case 21: more than one spectral lines, no Rp dependance, all lines are affected by the same wind
    - case 22: more than one spectral lines, no Rp dependance, all lines have same FWHM
    - case 23: more than one spectral lines, no Rp dependance, all lines are affected by the same wind and have same FWHM
    - case 24: more than one spectral lines, no Rp dependance, no winds
    - case 25: more than one spectral lines, no Rp dependance, no winds, all lines have same FWHM

                free_Rp     free_winds  shared_winds    shared_FWHM
    - case  0:  True        True        False           False               DEFAULT for single line
    - case  1:  True        False       False           False
    - case  2:  False       True        False           False
    - case  3:  False       False       False           False


    - case 10:  True        True        False           False               DEFAULT for multiple lines
    - case 11:  True        True        True            False
    - case 12:  True        True        False           True
    - case 13:  True        True        True            True
    - case 14:  True        False       False           False
    - case 15:  True        False       False           True

    - case 20:  False       True        False           False
    - case 21:  False       True        True            False
    - case 22:  False       True        False           True
    - case 23:  False       True        True            True
    - case 24:  False       False       False           False
    - case 25:  False       False       False           True
    """

    model_case = 10

    fit_pams = lines_dict['fit_parameters']
    # Added compativlity to "wrong" keys

    clv_rm_correction = lines_dict.get('clv_rm_correction', True)

    free_Rp = fit_pams.get('free_Rp', True) \
        and fit_pams.get('free_planet_radius', True) \
        and clv_rm_correction
    free_winds = fit_pams.get('free_winds', True) \
        and fit_pams.get('free_offset', True)
    shared_winds = fit_pams.get('shared_winds', False) \
        or fit_pams.get('shared_offset', False)
    shared_FWHM = fit_pams.get('shared_FWHM', False) \
        or fit_pams.get('shared_fwhm', False)

    prior_dict = fit_pams.get('priors', {}) \
        or fit_pams.get('priors', {})

    allow_emission = fit_pams.get('allow_emission', False)

    if len(lines_dict['lines']) < 2:
        if free_Rp is True and free_winds is True:
            model_case = 0
        if free_Rp is True and free_winds is False:
            model_case = 1
        if free_Rp is False and free_winds is True:
            model_case = 2
        if free_Rp is False and free_winds is False:
            model_case = 3
    else:
        if free_Rp is True:
            if free_winds is True:
                if shared_winds is False and shared_FWHM is False:
                    model_case = 10
                if shared_winds is True and shared_FWHM is False:
                    model_case = 11
                if shared_winds is False and shared_FWHM is True:
                    model_case = 12
                if shared_winds is True and shared_FWHM is True:
                    model_case = 13
            else:
                if shared_winds is False and shared_FWHM is False:
                    model_case = 14
                if shared_winds is False and shared_FWHM is True:
                    model_case = 15
        else:
            if free_winds is True:
                if shared_winds is False and shared_FWHM is False:
                    model_case = 20
                if shared_winds is True and shared_FWHM is False:
                    model_case = 21
                if shared_winds is False and shared_FWHM is True:
                    model_case = 22
                if shared_winds is True and shared_FWHM is True:
                    model_case = 23
            else:
                if shared_winds is False and shared_FWHM is False:
                    model_case = 24
                if shared_winds is False and shared_FWHM is True:
                    model_case = 25

    jitter_flag = fit_pams.get('jitter', True)

    print()
    print('      free_Rp:         (default: True)  ', free_Rp)
    print('      free_winds:      (default: True)  ', free_winds)
    print('      shared_winds:    (default: False) ', shared_winds)
    print('      shared_FWHM:     (default: False) ', shared_FWHM)
    print('      jitter:          (default: True)  ', jitter_flag)
    print('      # lines:        ', len(lines_dict['lines']))
    print('      model_case:     ', model_case)

    """ parameters list:
        to be updated


    pams_dict = {}  # dictionary containing the index of a given parameter
    pams_list = []  # list with the parameter names ordered according to their index
    boundaries = np.empty([0, 2])  # boundaries for MCMC / nested sampling
    theta_start = np.empty(0)  # starting point for MCMC
    lines_center = np.empty(0)  # laboratory wavelength of spectral lines
    pam_index = 0  # keep track of the number of variables

    for line_key, line_val in lines_dict['lines'].items():
        pam_name = line_key + '_contrast'
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        boundaries = np.append(boundaries, [[0.00, 1.00]], axis=0)
        theta_start = np.append(theta_start, 0.010)
        pam_index += 1

        lines_center = np.append(lines_center, line_val)

        # skip the inclusion of FWHM as a free parameter for each line
            if the shared FWHM is selected
        #

        if model_case in [0, 1, 2, 3, 10, 11, 14, 20, 21, 24]:
            # if not lines_dict['fit_parameters']['shared_fwhm']:
            pam_name = line_key + '_fwhm'
            pams_dict[pam_name] = pam_index
            pams_list.append(pam_name)
            boundaries = np.append(boundaries, [[0.00, 150.00]], axis=0)
            theta_start = np.append(theta_start, 5.0)
            pam_index += 1

        # if lines_dict['fit_parameters']['fixed_separation']: continue
        # if not lines_dict['fit_parameters']['lines_shift']: continue

        if model_case in [0, 2, 10, 12, 20, 22]:
            pam_name = line_key + '_winds'
            pams_dict[pam_name] = pam_index
            pams_list.append(pam_name)
            boundaries = np.append(boundaries, [[-5.00, 5.00]], axis=0)
            theta_start = np.append(theta_start, 0.00)
            pam_index += 1

    if model_case in [12, 13, 15, 22, 23, 25]:
        # if lines_dict['fit_parameters']['shared_fwhm']:
        pam_name = 'shared_fwhm'
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        boundaries = np.append(boundaries, [[0.000, 150.00]], axis=0)
        theta_start = np.append(theta_start, 5.000)
        pam_index += 1

    if model_case in [11, 13, 21, 23]:
        # if lines_dict['fit_parameters']['fixed_separation'] and lines_dict['fit_parameters']['lines_shift']:
        pam_name = 'shared_winds'
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        boundaries = np.append(boundaries, [[-5.0, 5.0]], axis=0)
        theta_start = np.append(theta_start, 0.000)
        pam_index += 1

    if model_case in [0, 1, 10, 11, 12, 13, 14, 15]:
        pams_dict['rp_factor'] = pam_index
        pams_list.append('rp_factor')
        boundaries = np.append(boundaries, [[0.5, 2.0]], axis=0)
        theta_start = np.append(theta_start, 1.0)
        pam_index += 1

    pams_dict['K_planet'] = pam_index
    pams_list.append('K_planet')
    boundaries = np.append(boundaries,
                            [[-300., planet_dict['RV_semiamplitude']
                                [0]+ 300.]],
                            axis=0)
    theta_start = np.append(
        theta_start, planet_dict['RV_semiamplitude'][0])
    pam_index += 1

    pam_name = 'jitter'
    pams_dict[pam_name] = pam_index
    pams_list.append(pam_name)
    boundaries = np.append(boundaries, [[10**(-12), 0.01]], axis=0)
    theta_start = np.append(theta_start, 10**(-11))
    pam_index += 1

    for ii in range(0, pam_index):
        print(pams_list[ii], '    ', boundaries[ii, :],
                '    ', theta_start[ii])

    ndim = pam_index

    """

    for night in night_dict:

        print()
        print("transmission_mcmc               Night: {0:s}".format(night))

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        observational_pams = load_from_cpickle(
            'observational_pams', config_in['output'], night)

        # Moved here to check wheter PCA or master-out have been employed
        preparation_input = load_from_cpickle(
            'transmission_preparation', config_in['output'], night)

        if preparation_input.get('pca_output', False):
            if pca_iteration >= 0:
                it_string = str(pca_iteration).zfill(2)
            else:
                it_string = str(preparation_input.get('ref_iteration', 0)).zfill(2)
            preparation = preparation_input[it_string]
        else:
            preparation = preparation_input
            it_string = ''

        try:
            mcmc_data = load_from_cpickle(subroutine_name + '_data', config_in['output'], night, lines_label, it_string)

            clv_rm_radius = mcmc_data['clv_rm_radius']
            clv_rm_grid = mcmc_data['clv_rm_grid']
            transmission_spec = mcmc_data['transmission_spec']
            transmission_spec_err = mcmc_data['transmission_spec_err']
            wave_array = mcmc_data['wave_array']
            time_array = mcmc_data['time_array']
            planet_RVsinusoid = mcmc_data['planet_RVsinusoid']
            jitter_index = mcmc_data['jitter_index']
            n_jitter = mcmc_data['n_jitter']

            print("   Loading MCMC data array for lines {0:s}, night: {1:s}".format(
                lines_label, night))

        except FileNotFoundError:

            print("   Computing MCMC data array for lines {0:s}, night: {1:s}".format(
                lines_label, night))

            calib_data = load_from_cpickle(
                'calibration_fibA', config_in['output'], night)
            input_data = retrieve_observations(
                config_in['output'], night, lists['observations'])

            if clv_rm_correction:
                try:
                    clv_rm_models = load_from_cpickle(
                        'clv_rm_models', config_in['output'], night, lines_label)
                except (FileNotFoundError, IOError):
                    clv_rm_models = load_from_cpickle(
                        'clv_rm_models', config_in['output'], night)
            else:
                # workaround if CLV correction is not available
                clv_rm_models = {'common': {}}
                clv_rm_models['common']['n_radius_grid'] = 3
                clv_rm_models['common']['radius_grid'] = np.asarray(
                    [0.5, 1.0, 1.5])

            processed = {
                'subroutine': subroutine_name,
            }

            """
            we use the first transit_full observation to define the boolean
            selection array to be used for all the observations.
            In this way we make sure that all the wavelength/spectrum arrays have the
            same dimension
            """
            obs_reference = lists['transit_full'][0]
            wave_SRF, step_SRF = shift_wavelength(input_data[obs_reference]['wave'],
                                                  input_data[obs_reference]['step'],
                                                  observational_pams[obs_reference]['rv_shift_ORF2SRF'])

            print('WAVE SHAPE:', np.shape(wave_SRF))
            print('STEP SHAPE:', np.shape(step_SRF))

            processed['common'] = {
                'range': lines_dict['fit_parameters']['range'],
                'reference_wave': wave_SRF,
                'reference_step': step_SRF
            }

            """
            Identification of the orders including the data points of interest
            """

            #identify_order = (wave_SRF > processed['common']['range'][0]) \
                #& (wave_SRF < processed['common']['range'][1])
            try:
                identify_order = load_from_cpickle(
                    'global_common_selection',
                    config_in['output'],
                    lines_label
                )
            except FileNotFoundError:
                identify_order = (wave_SRF > processed['common']['range'][0]) \
                              & (wave_SRF < processed['common']['range'][1])

            order_selection = (np.sum(identify_order, axis=1) > 0)
            order_list = np.arange(0, observational_pams['n_orders'], dtype=np.int16)[order_selection]
            n_orders = len(order_list)

            processed['common']['selection'] = identify_order
            processed['common']['order_selection'] = order_selection
            processed['common']['order_list'] = order_list
            processed['common']['n_orders'] = n_orders
            processed['common']['size'] = np.sum(identify_order)
            print("GLOBAL COMMON SIZE:", processed['common']['size'])

            # SAVE GLOBAL COMMON SELECTION (first night only)
            if night == list(night_dict.keys())[0]:
                save_to_cpickle(
                    'global_common_selection',
                    processed['common']['selection'],
                    config_in['output'],
                    lines_label
                )

            print('COMMON')
            print(np.shape(processed['common']['selection']))
            print(processed['common']['order_selection'])
            print(processed['common']['order_list'])
            print(processed['common']['size'])

            processed['common_extended'] = {
                'range': lines_dict['range'],
                'reference_wave': wave_SRF,
                'reference_step': step_SRF
            }

            identify_order = (wave_SRF > processed['common_extended']['range'][0]) \
                & (wave_SRF < processed['common_extended']['range'][1])
            order_selection = (np.sum(identify_order, axis=1) > 0)
            order_list = np.arange(0, observational_pams['n_orders'], dtype=np.int16)[order_selection]
            n_orders = len(order_list)

            processed['common_extended']['selection'] = identify_order
            processed['common_extended']['order_selection'] = order_selection
            processed['common_extended']['order_list'] = order_list
            processed['common_extended']['n_orders'] = n_orders
            processed['common_extended']['size'] = np.sum(identify_order)

            #print('COMMON_EXTENDED')
            #print(np.shape(processed['common_extended']['selection']))
            #print(processed['common_extended']['order_selection'])
            #print(processed['common_extended']['order_list'])
            #print(processed['common_extended']['size'])
            #print()

            for obs in lists['observations']:

                """ we start from the e2ds file, after correction for blaze and
                    division by the master-out

                Observation data:
                    wave: input_data[obs]['wave']
                    step: input_data[obs]['step']
                    flux: preparation[obs]['deblazed']
                    ferr: preparation[obs]['deblazed_err']
                """

                """ First step: we rebin the spectra in the Stellar Reference Frame,
                    with the step size decided by the user specifically for the fit
                """
                processed[obs] = {}
                processed[obs]['wave_SRF'], processed[obs]['step_SRF'] = shift_wavelength(
                    input_data[obs]['wave'],
                    input_data[obs]['step'],
                    observational_pams[obs]['rv_shift_ORF2SRF'])


                """ Continuum normalization:
                    3) Polynomial fit, everything is hard coded now but personalized
                    options can be implemented easily in the yaml file
                """

                processed[obs]['continuum'] = processed['common']['reference_wave'] * 0.
                processed[obs]['normalized'] = processed['common']['reference_wave'] * 0.
                processed[obs]['normalized_err'] = processed['common']['reference_wave'] * 0.
                """ We perform the selection only on the order
                    where we actually have data for the MCMC
                """




            if norm_pams['normalize_transmission'] and norm_pams['normalization_model'] == 'polynomial':

                """ Continuum normalization preparatory steps:
                    1) exclusion of regions with planetary lines
                    2) exclusion of regions with stellar lines
                    3) Polynomial fit of selected regions
                    Boolean array initialized to all True values, fit is
                    performed on the extended region and then applied to the fit subset
                """
                processed['common_extended']['line_exclusion'] = processed['common_extended']['selection'].copy()

                """ Continuum normalization:
                    1) exclusion of regions with planetary lines, taking into
                    account the planetary RV semi-amplitude
                """
                for line_key, line_val in lines_dict['lines'].items():
                    line_extension = 1.2 * \
                        planet_dict['RV_semiamplitude'][0] * \
                        line_val / speed_of_light_km
                    processed['common_extended']['line_exclusion'] = processed['common_extended']['line_exclusion'] \
                        & (np.abs(processed['common_extended']['reference_wave']-line_val) > line_extension)

                """ Continuum normalization:
                    2) exclusion of regions with planetary lines, taking into
                    account the planetary RV semi-amplitude
                """

                try:

                    for order in processed['common_extended']['order_list']:
                        print('ORDER', order)
                        order_sel = processed['common_extended']['selection'][order, :]
                        print('selection', np.shape(processed['common_extended']['selection']))
                        print('ORDER_SEL', np.shape(order_sel))
                        stellar_spectrum_rebinned = rebin_1d_to_1d(clv_rm_models['common']['wave'],
                                                                clv_rm_models['common']['step'],
                                                                clv_rm_models['common']['norm_convolved'],
                                                                processed['common_extended']['reference_wave'][order, order_sel],
                                                                processed['common_extended']['reference_step'][order, order_sel])

                        stellar_spectrum_derivative = first_derivative(
                            processed['common_extended']['reference_wave'][order, order_sel], stellar_spectrum_rebinned)

                        processed['common_extended']['line_exclusion'][order, order_sel] = processed['common_extended']['line_exclusion'][order, order_sel] & (
                            np.abs(stellar_spectrum_derivative) < 0.0005)
                except KeyError:
                    print(
                        "No stellar synthetic spectrum from CLV models, some stellar lines may be included transmission normalization  ")

                for obs in lists['observations']:
                    for order in processed['common']['order_list']:

                        selection = processed['common_extended']['line_exclusion'][order, :] & (
                            preparation[obs]['deblazed'][order, :]
                            > np.std(preparation[obs]['deblazed'][order, :]))

                        if np.sum(selection) < 100:
                            selection = (preparation[obs]['deblazed'][order, :] >
                                        np.std(preparation[obs]['deblazed'][order, :]))

                        processed[obs]['norm_coeff_' + repr(order)] = \
                            np.polynomial.chebyshev.chebfit(processed[obs]['wave_SRF'][order, selection],
                                                            preparation[obs]['deblazed'][order, selection],
                                                            2)
                        processed[obs]['continuum'][order, selection] = np.polynomial.chebyshev.chebval(
                            processed[obs]['wave_SRF'][order, selection], processed[obs]['norm_coeff_' + repr(order)])
                        processed[obs]['normalized'][order, selection] = preparation[obs]['deblazed'][order, selection] / \
                            processed[obs]['continuum'][order, selection]
                        processed[obs]['normalized_err'][order, selection] = preparation[obs]['deblazed_err'][order, selection] / \
                            processed[obs]['continuum'][order, selection]


            elif norm_pams['normalize_transmission']  and (
                norm_pams['normalization_model'] == 'savgol'
                or norm_pams['normalization_model'] == 'savitzky-golay'):

                print(' ', obs, ' normalization using Savitzky-Golay filter')



                for obs in lists['observations']:

                    processed[obs]['continuum'] = np.pnes_like(preparation[obs]['deblazed'])
                    for order in processed['common']['order_list']:

                        processed[obs]['continuum'][order,:] = savgol_filter(preparation[obs]['deblazed'][order,:],
                                                                window_length=norm_pams['window_length'],
                                                                polyorder=norm_pams['polyorder'],
                                                                mode=norm_pams['mode'],
                                                                cval=norm_pams['cval'])

                    processed[obs]['normalized'] = preparation[obs]['deblazed'] /  processed[obs]['continuum']
                    processed[obs]['normalized_err'] = preparation[obs]['deblazed_err'] / processed[obs]['continuum']


                ######
                print('ciao')





            processed['common']['n_obs'] = len(lists['transit_full'])
            processed['common']['n_radius_grid'] = clv_rm_models['common']['n_radius_grid']
            processed['common']['radius_grid'] = clv_rm_models['common']['radius_grid']
            clv_rm_radius = clv_rm_models['common']['radius_grid']

            """ We are moving the values of interest from dictionaries to arrays
                in order to  speed up the MCMC
                1) clv_rm_grid: array with all the CLV models, as a function of the
                    radius of the planet
                2) time_from_transit: BJD_TDB - T0
                3) planet_RVsinusoid: Fractional RV of the planet (K=1) - from a array
            """
            clv_rm_grid = np.ones([processed['common']['n_radius_grid'],
                                   processed['common']['n_obs'],
                                   processed['common']['size']],
                                  dtype=np.double)

            time_from_transit = np.empty(
                processed['common']['n_obs'], dtype=np.double)

            wave_array = np.empty([processed['common']['n_obs'],
                                   processed['common']['size']],
                                  dtype=np.double)
            time_array = np.empty([processed['common']['n_obs'],
                                   processed['common']['size']],
                                  dtype=np.double)
            transmission_spec = np.empty([processed['common']['n_obs'],
                                          processed['common']['size']],
                                         dtype=np.double)
            transmission_spec_err = np.empty([processed['common']['n_obs'],
                                              processed['common']['size']],
                                             dtype=np.double)

            for i_obs, obs in enumerate(lists['transit_full']):

                time_from_transit[i_obs] = observational_pams[obs]['BJD'] - \
                    observational_pams['time_of_transit']
                time_array[i_obs, :] = time_from_transit[i_obs]

                # planet_RVsinusoid[i_obs] = np.sin(2*np.pi /  planet_dict['period'][0] * time_from_transit[i_obs])
                wave_array[i_obs, :] = processed[obs]['wave_SRF'][processed['common']['selection']].flatten()
                transmission_spec[i_obs, :] = processed[obs]['normalized'][processed['common']['selection']].flatten()
                transmission_spec_err[i_obs,
                                      :] = processed[obs]['normalized_err'][processed['common']['selection']].flatten()

                if clv_rm_correction is False:
                    continue

                for i_r in range(0, processed['common']['n_radius_grid']):

                    clv_rm_temp = processed['common_extended']['reference_wave'] * 0.
                    for order in processed['common']['order_list']:

                        """ CLV Synthetic models are in the Stellar Reference system,
                        so no shift is required """
                        clv_rm_temp[order, :] = rebin_1d_to_1d(
                            clv_rm_models['common']['wave'],
                            clv_rm_models['common']['step'],
                            clv_rm_models[obs]['clv_rm_model_convolved_normalized'][i_r, :],
                            processed[obs]['wave_SRF'][order, :],
                            processed[obs]['step_SRF'][order, :],
                            preserve_flux=False)

                    clv_rm_grid[i_r, i_obs, :] = clv_rm_temp[processed['common']['selection']].flatten()

                    # preserve_flux should be True or False?
                    # False if the spectra are already normalized

            remove_outliers = (np.abs(transmission_spec - 1.) > 0.5)
            transmission_spec[remove_outliers] = 1.0
            transmission_spec_err[remove_outliers] = 1.0

            planet_RVsinusoid = np.sin(
                2*np.pi / planet_dict['period'][0] * time_array)

            if jitter_flag:
                jitter_index = []
                n_jitter = 1
            else:
                jitter_index = None
                n_jitter = 0

            mcmc_data = {
                'observations': lists['transit_full'],
                'common_wave': processed['common']['reference_wave'],
                'common_step': processed['common']['reference_step'],
                'clv_rm_grid': clv_rm_grid,
                'transmission_spec': transmission_spec,
                'transmission_spec_err': transmission_spec_err,
                'wave_array': wave_array,
                'time_array': time_array,
                'planet_RVsinusoid': planet_RVsinusoid,
                'clv_rm_radius': clv_rm_models['common']['radius_grid'],
                'n_obs': len(lists['transit_full']),
                'n_radius_grid': clv_rm_models['common']['n_radius_grid'],
                'jitter_index': jitter_index,
                'n_jitter': n_jitter
            }

            save_to_cpickle(subroutine_name + '_data', mcmc_data,
                            config_in['output'], night, lines_label, it_string)

            # Forcing memory deallocation
            clv_rm_models = None
            mcmc_data = None

        print()
        print("transmission_mcmc           ")

        try:
            results_dict = load_from_cpickle(subroutine_name+'_'+sampler_name+'_results',
                                             config_in['output'], night, lines_label, it_string)
            print("   Transmission MCMC analysis for lines {0:s}, night: {1:s}  already performed".format(
                lines_label, night))

            pams_dict = results_dict['pams_dict']
            chain_med = results_dict['chain_med']
            boundaries = results_dict['boundaries']
            start_average = np.average(results_dict['point_start'], axis=0)
            ndim = results_dict['ndim']
            # TODO improve output
            print('   *** sampler output ')

            for key, val in pams_dict.items():
                print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}]) (start: {7:9f})'.format(key, val,
                                                                                                                        chain_med[val, 0],
                                                                                                                        chain_med[val, 2],
                                                                                                                        chain_med[val, 1],
                                                                                                                        boundaries[val, 0],
                                                                                                                        boundaries[val, 1],
                                                                                                                        start_average[val])
                      )

            continue

            # R(h) = np.sqrt(1+h/delta)

        except FileNotFoundError:
            print()

        # getting fit parameters
        lines_center, pams_dict, pams_list, boundaries, theta_start = define_theta_array(
            model_case, lines_dict, planet_dict, n_jitter,allow_emission=allow_emission)
        ndim = len(theta_start)
        ngen = sampler_pams.get('n_gen', 64000)
        nwalkers_mult = sampler_pams.get('n_walkers_mult', 2)
        nwalkers = sampler_pams.get('n_walkers', nwalkers_mult * ndim)
        nthin = sampler_pams.get('n_thin', 50)
        nsteps = sampler_pams.get('n_steps', 20000)
        nburnin = sampler_pams.get('n_burnin', 5000)
        ncpus = sampler_pams.get('n_cpus', 4)
        ndata = np.size(wave_array)

        if pams_dict.get('rp_factor', False):
            pam_id = pams_dict['rp_factor']
            boundaries[pam_id, :] = [clv_rm_radius[0], clv_rm_radius[-1]]

        print()
        print('      PyDE + emcee parameters')
        print('         n_dim:                        {0:9.0f}'.format(ndim))
        print(
            '         n_walkers:  (default: 2*ndim) {0:9.0f}'.format(nwalkers))
        print('         n_gen:      (default: 64000)  {0:9.0f}'.format(ngen))
        print('         n_steps:    (default: 20000)  {0:9.0f}'.format(nsteps))
        print(
            '         n_burnin:   (default: 10000)  {0:9.0f}'.format(nburnin))
        print('         n_thin:     (default: 50)     {0:9.0f}'.format(nthin))

        population, sampler_chain, sampler_lnprobability, point_start = emcee_lines_fit_functions(
            model_case,
            wave_array,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            prior_dict,
            theta_start, boundaries, ndim, nwalkers, ngen, nsteps, nthin, ncpus)

        flat_chain, flat_lnprob, chain_med, chain_MAP, lnprob_med, lnprob_MAP = \
            emcee_flatten_median(population, sampler_chain,
                                 sampler_lnprobability,  nburnin, nthin, nwalkers)
        emcee_compute_BIC_AIC(lnprob_med, lnprob_MAP, ndata, ndim)

        med_lines_model, med_clv_model, med_lines_array, med_planet_K, med_planet_R, med_jitter = \
            return_model(model_case,
                         chain_med[:, 0],
                         wave_array,
                         clv_rm_radius,
                         clv_rm_grid,
                         planet_RVsinusoid,
                         lines_center,
                         jitter_index)

        map_lines_model, map_clv_model, map_lines_array, map_planet_K, map_planet_R, map_jitter = \
            return_model(model_case,
                         chain_MAP,
                         wave_array,
                         clv_rm_radius,
                         clv_rm_grid,
                         planet_RVsinusoid,
                         lines_center,
                         jitter_index)

        results_dict = {
            'sampler_name': sampler_name,
            'ndim': ndim,
            'nwalkers': nwalkers,
            'nthin': nthin,
            'nsteps': nsteps,
            'nburnin': nburnin,
            'ndata': ndata,
            'pams_dict': pams_dict,
            'population': population,
            'sampler_chain': sampler_chain,
            'sampler_lnprobability': sampler_lnprobability,
            'theta_start': theta_start,
            'boundaries': boundaries,
            'flat_chain': flat_chain,
            'flat_lnprob': flat_lnprob,
            'chain_med': chain_med,
            'chain_MAP': chain_MAP,
            'lnprob_med': lnprob_med,
            'lnprob_MAP': lnprob_MAP,
            'lines_center': lines_center,
            'point_start': point_start,
            'theta_start': theta_start,
        }

        results_dict['results'] = {
            'lines_model': med_lines_model,
            'clv_model': med_clv_model,
            'lines_array': med_lines_array,
            'planet_K': med_planet_K,
            'planet_R': med_planet_R,
            'jitter': med_jitter
        }

        results_dict['results_MAP'] = {
            'lines_model': map_lines_model,
            'clv_model': map_clv_model,
            'lines_array': map_lines_array,
            'planet_K': map_planet_K,
            'planet_R': map_planet_R,
            'jitter': map_jitter
        }

        results_dict['results']['observational_pams'] = {}
        results_dict['results_MAP']['observational_pams'] = {}

        for obs in lists['observations']:

            results_dict['results']['observational_pams'][obs] = {}
            results_dict['results_MAP']['observational_pams'][obs] = {}

            """ RV shift from the observer RF to the planet RF
                STRONG ASSUMPTIONS:
                    - there is only the transiting planet in the system
                    - the planet has null eccentricity
                    - linear approximation or the orbit near the transit event

                Computation is performed by moving to the Solar Barycenter, than to the Stellar System Barycenter
                and finally onto the planet
            """
            results_dict['results']['observational_pams'][obs]['rv_shift_ORF2PRF'] = \
                observational_pams[obs]['BERV'] \
                - observational_pams['RV_star']['RV_systemic'] \
                - results_dict['results']['planet_K'] \
                * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                / planet_dict['period'][0] * 2 * np.pi

            results_dict['results_MAP']['observational_pams'][obs]['rv_shift_ORF2PRF'] = \
                observational_pams[obs]['BERV'] \
                - observational_pams['RV_star']['RV_systemic'] \
                - results_dict['results_MAP']['planet_K'] \
                * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                / planet_dict['period'][0] * 2 * np.pi

            """ RV shift from Stellar Rest Frame to Planetary Rest Frame
                We have to take into account the RV of star relatively to the Barycenter
            """
            results_dict['results']['observational_pams'][obs]['rv_shift_SRF2PRF'] = \
                + observational_pams[obs]['RV_bjdshift'] \
                - results_dict['results']['planet_K'] \
                * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                / planet_dict['period'][0] * 2 * np.pi

            results_dict['results_MAP']['observational_pams'][obs]['rv_shift_SRF2PRF'] = \
                + observational_pams[obs]['RV_bjdshift'] \
                - results_dict['results_MAP']['planet_K'] \
                * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                / planet_dict['period'][0] * 2 * np.pi

        save_to_cpickle(subroutine_name+'_'+sampler_name+'_results',
                        results_dict, config_in['output'], night, lines_label, it_string)

        # TODO improve output
        print('   *** sampler output ')

        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}])'.format(key, val,
                                                                                                    chain_med[val, 0],
                                                                                                    chain_med[val, 2],
                                                                                                    chain_med[val, 1],
                                                                                                    boundaries[val, 0], boundaries[val, 1])
                  )
        # print('   *** physical output')
        #
        #        results_dict['results'] = {
        #    'lines_model': med_lines_model,
        #    'clv_model': med_clv_model,
        #    'lines_array': med_lines_array,
        #    'planet_K': med_planet_K,
        #    'planet_R': med_planet_R,
        #    'jitter': med_jitter
        # }

    """ Analysis of the entire dataset """

    print()
    try:
        all_mcmc_data = load_from_cpickle(
            subroutine_name+'_data', config_in['output'], night='', lines=lines_label, it_string=it_string)

        all_clv_rm_radius = all_mcmc_data['clv_rm_radius']
        all_clv_rm_grid = all_mcmc_data['clv_rm_grid']
        all_transmission_spec = all_mcmc_data['transmission_spec']
        all_transmission_spec_err = all_mcmc_data['transmission_spec_err']
        all_wave_array = all_mcmc_data['wave_array']
        all_time_array = all_mcmc_data['time_array']
        all_planet_RVsinusoid = all_mcmc_data['planet_RVsinusoid']
        all_observations = all_mcmc_data['observations']
        all_n_obs = all_mcmc_data['n_obs']
        all_n_radius_grid = all_mcmc_data['n_radius_grid']
        all_jitter_index = all_mcmc_data['jitter_index']
        n_jitter = all_mcmc_data['n_jitter']

    except:

        n_jitter = 0

        for night in night_dict:

            mcmc_data = load_from_cpickle(subroutine_name+'_data',
                                          config_in['output'], night, lines_label, it_string=it_string)

            try:
                # Building the arrays for the full analysis

                all_clv_rm_grid = np.concatenate(
                    (all_clv_rm_grid, mcmc_data['clv_rm_grid']), axis=1)
                all_transmission_spec = np.concatenate(
                    (all_transmission_spec, mcmc_data['transmission_spec']))
                all_transmission_spec_err = np.concatenate(
                    (all_transmission_spec_err, mcmc_data['transmission_spec_err']))
                all_wave_array = np.concatenate(
                    (all_wave_array, mcmc_data['wave_array']))
                all_time_array = np.concatenate(
                    (all_time_array, mcmc_data['time_array']))
                all_planet_RVsinusoid = np.concatenate(
                    (all_planet_RVsinusoid, mcmc_data['planet_RVsinusoid']))
                all_observations = np.concatenate(
                    (all_observations, mcmc_data['observations']))
                all_n_obs += mcmc_data['n_obs']

                if jitter_flag:
                    all_jitter_index = np.concatenate(
                        (all_jitter_index, n_jitter*np.ones(np.shape(mcmc_data['wave_array']), dtype=np.int16)))
                    n_jitter += 1

            except NameError:
                """ This error is expected when retrieving the data of the first night"""
                all_clv_rm_radius = mcmc_data['clv_rm_radius']
                all_clv_rm_grid = mcmc_data['clv_rm_grid']
                all_transmission_spec = mcmc_data['transmission_spec']
                all_transmission_spec_err = mcmc_data['transmission_spec_err']
                all_wave_array = mcmc_data['wave_array']
                all_time_array = mcmc_data['time_array']
                all_planet_RVsinusoid = mcmc_data['planet_RVsinusoid']
                all_observations = mcmc_data['observations']
                all_n_obs = mcmc_data['n_obs']
                all_n_radius_grid = mcmc_data['n_radius_grid']
                if jitter_flag:
                    all_jitter_index = n_jitter * \
                        np.ones(
                            np.shape(mcmc_data['wave_array']), dtype=np.int16)
                    n_jitter += 1
                else:
                    all_jitter_index = None

        all_mcmc_data = {
            'observations': all_observations,
            'clv_rm_grid': all_clv_rm_grid,
            'transmission_spec': all_transmission_spec,
            'transmission_spec_err': all_transmission_spec_err,
            'wave_array': all_wave_array,
            'time_array': all_time_array,
            'planet_RVsinusoid': all_planet_RVsinusoid,
            'clv_rm_radius': all_clv_rm_radius,
            'n_obs': all_n_obs,
            'n_radius_grid': all_n_radius_grid,
            'jitter_index': all_jitter_index,
            'n_jitter': n_jitter
        }

        save_to_cpickle(subroutine_name+'_data', all_mcmc_data,
                        config_in['output'], night='', lines=lines_label, it_string=it_string)

    try:
        results_dict = load_from_cpickle(subroutine_name + '_' + sampler_name+'_results',
                                         config_in['output'], night='', lines=lines_label, it_string=it_string)
        print("   Transmission MCMC analysis for lines {0:s} already performed ".format(
            lines_label))

        pams_dict = results_dict['pams_dict']
        chain_med = results_dict['chain_med']
        boundaries = results_dict['boundaries']
        ndim = results_dict['ndim']
        # TODO improve output
        print('   *** sampler output ')

        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}])'.format(key, val,
                                                                                                    chain_med[val, 0],
                                                                                                    chain_med[val, 2],
                                                                                                    chain_med[val, 1],
                                                                                                    boundaries[val, 0], boundaries[val, 1])
                  )

    except FileNotFoundError:

        lines_center, pams_dict, pams_list, boundaries, theta_start = define_theta_array(
            model_case, lines_dict, planet_dict, n_jitter,allow_emission=allow_emission)
        ndim = len(theta_start)

        if pams_dict.get('rp_factor', False):
            pam_id = pams_dict['rp_factor']
            boundaries[pam_id, :] = [clv_rm_radius[0], clv_rm_radius[-1]]

        ndata = np.size(all_wave_array)

        print()
        print('      PyDE + emcee parameters')
        print('         n_dim:                        {0:9.0f}'.format(ndim))
        print(
            '         n_walkers:  (default: 2*ndim) {0:9.0f}'.format(nwalkers))
        print('         n_gen:      (default: 64000)  {0:9.0f}'.format(ngen))
        print('         n_steps:    (default: 20000)  {0:9.0f}'.format(nsteps))
        print(
            '         n_burnin:   (default: 10000)  {0:9.0f}'.format(nburnin))
        print('         n_thin:     (default: 50)     {0:9.0f}'.format(nthin))

        population, sampler_chain, sampler_lnprobability, point_start = emcee_lines_fit_functions(
            model_case,
            all_wave_array,
            all_transmission_spec,
            all_transmission_spec_err,
            all_clv_rm_radius,
            all_clv_rm_grid,
            all_planet_RVsinusoid,
            lines_center,
            all_jitter_index,
            prior_dict,
            theta_start, boundaries, ndim, nwalkers, ngen, nsteps, nthin, ncpus)

        flat_chain, flat_lnprob, chain_med, chain_MAP, lnprob_med, lnprob_MAP = \
            emcee_flatten_median(population, sampler_chain,
                                 sampler_lnprobability,  nburnin, nthin, nwalkers)
        emcee_compute_BIC_AIC(lnprob_med, lnprob_MAP, ndata, ndim)

        med_lines_model, med_clv_model, med_lines_array, med_planet_K, med_planet_R, med_jitter = \
            return_model(model_case,
                         chain_med[:, 0],
                         all_wave_array,
                         all_clv_rm_radius,
                         all_clv_rm_grid,
                         all_planet_RVsinusoid,
                         lines_center,
                         all_jitter_index)

        map_lines_model, map_clv_model, map_lines_array, map_planet_K, map_planet_R, map_jitter = \
            return_model(model_case,
                         chain_MAP,
                         all_wave_array,
                         all_clv_rm_radius,
                         all_clv_rm_grid,
                         all_planet_RVsinusoid,
                         lines_center,
                         all_jitter_index)

        results_dict = {
            'sampler_name': sampler_name,
            'ndim': ndim,
            'nwalkers': nwalkers,
            'nthin': nthin,
            'nsteps': nsteps,
            'nburnin': nburnin,
            'ndata': ndata,
            'pams_dict': pams_dict,
            'population': population,
            'sampler_chain': sampler_chain,
            'sampler_lnprobability': sampler_lnprobability,
            'theta_start': theta_start,
            'boundaries': boundaries,
            'flat_chain': flat_chain,
            'flat_lnprob': flat_lnprob,
            'chain_med': chain_med,
            'chain_MAP': chain_MAP,
            'lnprob_med': lnprob_med,
            'lnprob_MAP': lnprob_MAP,
            'lines_center': lines_center,
            'point_start': point_start,
            'theta_start': theta_start,
        }

        results_dict['results'] = {
            'lines_model': med_lines_model,
            'clv_model': med_clv_model,
            'lines_array': med_lines_array,
            'planet_K': med_planet_K,
            'planet_R': med_planet_R,
            'jitter': med_jitter
        }

        results_dict['results_MAP'] = {
            'lines_model': map_lines_model,
            'clv_model': map_clv_model,
            'lines_array': map_lines_array,
            'planet_K': map_planet_K,
            'planet_R': map_planet_R,
            'jitter': map_jitter
        }

        results_dict['results']['observational_pams'] = {}
        results_dict['results_MAP']['observational_pams'] = {}

        for night in night_dict:

            lists = load_from_cpickle('lists', config_in['output'], night)
            observational_pams = load_from_cpickle(
                'observational_pams', config_in['output'], night)

            """ No differentiation by night """
            for obs in lists['observations']:

                results_dict['results']['observational_pams'][obs] = {}
                results_dict['results_MAP']['observational_pams'][obs] = {}

                """ RV shift from the observer RF to the planet RF
                    STRONG ASSUMPTIONS:
                        - there is only the transiting planet in the system
                        - the planet has null eccentricity
                        - linear approximation or the orbit near the transit event

                    Computation is performed by moving to the Solar Barycenter, than to the Stellar System Barycenter
                    and finally onto the planet
                """
                results_dict['results']['observational_pams'][obs]['rv_shift_ORF2PRF'] = \
                    observational_pams[obs]['BERV'] \
                    - observational_pams['RV_star']['RV_systemic'] \
                    - results_dict['results']['planet_K'] \
                    * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                    / planet_dict['period'][0] * 2 * np.pi

                results_dict['results_MAP']['observational_pams'][obs]['rv_shift_ORF2PRF'] = \
                    observational_pams[obs]['BERV'] \
                    - observational_pams['RV_star']['RV_systemic'] \
                    - results_dict['results_MAP']['planet_K'] \
                    * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                    / planet_dict['period'][0] * 2 * np.pi

                """ RV shift from Stellar Rest Frame to Planetary Rest Frame
                    We have to take into account the RV of star relatively to the Barycenter
                """
                results_dict['results']['observational_pams'][obs]['rv_shift_SRF2PRF'] = \
                    + observational_pams[obs]['RV_bjdshift'] \
                    - results_dict['results']['planet_K'] \
                    * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                    / planet_dict['period'][0] * 2 * np.pi

                results_dict['results_MAP']['observational_pams'][obs]['rv_shift_SRF2PRF'] = \
                    + observational_pams[obs]['RV_bjdshift'] \
                    - results_dict['results_MAP']['planet_K'] \
                    * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
                    / planet_dict['period'][0] * 2 * np.pi

        save_to_cpickle(subroutine_name + '_'+sampler_name+'_results',
                        results_dict, config_in['output'], night='', lines=lines_label, it_string=it_string)

        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}])'.format(key, val,
                                                                                                    chain_med[val, 0],
                                                                                                    chain_med[val, 2],
                                                                                                    chain_med[val, 1],
                                                                                                    boundaries[val, 0], boundaries[val, 1])
                  )
        print('MCMC completed')

    # Update planet parameters
    # deprecated
    # try:
    #    _ = load_from_cpickle(
    #        'observational', config_in['output'], night, lines_label)
    #    print("   Transmission MCMC results for lines {0:s} already store in observational array".format(
    #            lines_label))
    # except FileNotFoundError:
    #
    #    results_full = load_from_cpickle('transmission_mcmc_'+sampler_name+'_results',
    #                                     config_in['output'], night='', lines=lines_label)
    #
    #    for night in night_dict:
    #
    #        results_night = load_from_cpickle('transmission_mcmc_'+sampler_name+'_results',
    #                                     config_in['output'], night=night, lines=lines_label)

    #        lists = load_from_cpickle('lists', config_in['output'], night)
    #        observational_pams = load_from_cpickle(
    #            'observational_pams', config_in['output'], night)

    #        for obs in lists['observations']:
    #
    #            """ RV shift from the observer RF to the planet RF
    #                STRONG ASSUMPTIONS:
    #                    - there is only the transiting planet in the system
    #                    - the planet has null eccentricity
    #                    - linear approximation or the orbit near the transit event
    #
    #                Computation is performed by moving to the Solar Barycenter, than to the Stellar System Barycenter
    #                and finally onto the planet
    #            """
    #            observational_pams[obs]['rv_shift_ORF2PRF'] = \
    #                observational_pams[obs]['BERV'] \
    #                - observational_pams['RV_star']['RV_systemic'] \
    #                - results_full['results']['planet_K'] \
    #                * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
    #                / planet_dict['period'][0] * 2 * np.pi

    #            """ RV shift from Stellar Rest Frame to Planetary Rest Frame
    #                We have to take into account the RV of star relatively to the Barycenter
    #            """
    #            observational_pams[obs]['rv_shift_SRF2PRF'] = \
    #                + observational_pams[obs]['RV_bjdshift'] \
    #                - results_full['results']['planet_K'] \
    #                * (observational_pams[obs]['BJD'] - observational_pams['time_of_transit']) \
    #                / planet_dict['period'][0] * 2 * np.pi

    #        observational_pams['Rp_factor'] = results_full['results']['planet_R']
    #        observational_pams['lines_array'] = results_full['results']['lines_array']
    #        observational_pams['jitter'] = results_full['results']['jitter']

    #        save_to_cpickle('observational', observational_pams,
    #                        config_in['output'], night, lines_label)
