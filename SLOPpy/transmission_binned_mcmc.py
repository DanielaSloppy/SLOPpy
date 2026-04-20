from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.constants import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.math_functions import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.bayesian_emcee import *
from SLOPpy.subroutines.kepler_exo import compute_planet_RV
# from SLOPpy.subroutines.rebin_subroutines import *
from scipy.signal import savgol_filter

__all__ = ['compute_transmission_binned_mcmc','compute_transmission_binned_mcmc_iterative',
           'plot_transmission_binned_mcmc','plot_transmission_binned_mcmc_iterative']

subroutine_name = 'transmission_binned_mcmc'


def compute_transmission_binned_mcmc_iterative(config_in, lines_label):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0, pca_parameters.get('iterations',5)):
        compute_transmission_binned_mcmc(config_in, lines_label, reference='planetRF', pca_iteration=it)



def plot_transmission_binned_mcmc_iterative(config_in, lines_label, night_input='', reference='planetRF'):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0, pca_parameters.get('iterations',5)):
        plot_transmission_binned_mcmc(config_in, lines_label, night_input=night_input, reference=reference, pca_iteration=it)





def compute_transmission_binned_mcmc(config_in, lines_label, reference='planetRF', pca_iteration=-1):

    night_dict = from_config_get_nights(config_in)
    planet_dict = from_config_get_planet(config_in)
    star_dict = from_config_get_star(config_in)
    clv_rm_dict = from_config_get_clv_rm(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    """ selection of those parameters that are specific of the spectral line(s)
        under analysis
    """

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

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

    norm_dict = lines_dict.get('normalization', clv_rm_dict.get('normalization', {}))
    norm_pams={}
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
    norm_pams['normalize_rebinned'] = norm_dict.get('normalize_rebinned', False)

    # Added back-compatibility to old or "wrong" keys
    clv_rm_correction = lines_dict.get('clv_rm_correction', True)

    fit_pams = lines_dict['fit_parameters']
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
    pyde_flag = fit_pams.get('pyde', True)

    print()
    print('      free_Rp:         (default: True)  ', free_Rp)
    print('      free_winds:      (default: True) ', free_winds)
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

        ### Moved here to retrieve infomration about wheter PCA correction has been performed or not
        preparation_input = load_from_cpickle('transmission_preparation', config_in['output'], night)

        if preparation_input.get('pca_output', False):
            if pca_iteration >= 0:
                it_string = str(pca_iteration).zfill(2)
            else:
                it_string = str(preparation_input.get('ref_iteration', 0)).zfill(2)
            preparation = preparation_input[it_string]
        else:
            preparation = preparation_input
            it_string = ''

        """ This need to be checked only once, so it's ok to take the output of the last night
            and propagate it to the rest of subroutine
        """
        if len(it_string) > 0:
            pca_output = True
        else:
            pca_output = False

        try:
            mcmc_data = load_from_cpickle(subroutine_name + '_data', config_in['output'], night, lines_label, it_string)

            clv_rm_radius = mcmc_data['clv_rm_radius']
            clv_rm_grid = mcmc_data['clv_rm_grid']
            transmission_spec = mcmc_data['transmission_spec']
            transmission_spec_err = mcmc_data['transmission_spec_err']
            wave_meshgrid = mcmc_data['wave_meshgrid']
            time_meshgrid = mcmc_data['time_meshgrid']
            planet_RVsinusoid = mcmc_data['planet_RVsinusoid']
            jitter_index = mcmc_data['jitter_index']
            n_jitter = mcmc_data['n_jitter']

            print("   Loading MCMC data array for lines {0:s}, night: {1:s}".format(
                lines_label, night))

        except FileNotFoundError:

            print("  Computing MCMC data array for lines {0:s}, night: {1:s}".format(
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

            processed['common'] = {
                'range': lines_dict['fit_parameters']['range']
            }
            processed['common']['wave'] = np.arange(processed['common']['range'][0],
                                                    processed['common']['range'][1],
                                                    lines_dict['fit_parameters']['bin_step'],
                                                    dtype=np.double)
            processed['common']['size'] = len(processed['common']['wave'])
            processed['common']['step'] = np.ones(
                processed['common']['size'], dtype=np.double) * lines_dict['fit_parameters']['bin_step']

            processed['common_extended'] = {
                'range': lines_dict['range']
            }
            processed['common_extended']['wave'] = np.arange(processed['common_extended']['range'][0],
                                                    processed['common_extended']['range'][1],
                                                    lines_dict['fit_parameters']['bin_step'],
                                                    dtype=np.double)
            processed['common_extended']['size'] = len(
                processed['common_extended']['wave'])
            processed['common_extended']['step'] = np.ones(
                processed['common_extended']['size'], dtype=np.double) * lines_dict['fit_parameters']['bin_step']

            for obs in lists['observations']:

                """ we start from the e2ds file, after correction for blaze and
                    division by the master-out

                Observation data:
                    wave: input_data[obs]['wave']
                    step: input_data[obs]['step']
                    flux: preparation[obs]['ratio']
                    ferr: preparation[obs]['ratio']
                """

                """ First step: we rebin the spectra in the Stellar Reference Frame,
                    with the step size decided by the user specifically for the fit
                """

                if pca_output:
                    preserve_flux = False
                    blaze = np.ones_like(calib_data['blaze'])
                else:
                    preserve_flux = input_data[obs].get('absolute_flux', True)
                    blaze = calib_data['blaze']

                processed[obs] = {}
                processed[obs]['rebinned'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio'],
                                    blaze,
                                    processed['common']['wave'],
                                    processed['common']['step'],
                                    rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'],
                                    preserve_flux=preserve_flux)
                processed[obs]['rebinned_err'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio_err'],
                                    blaze,
                                    processed['common']['wave'],
                                    processed['common']['step'],
                                    rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'],
                                    preserve_flux=preserve_flux,
                                    is_error=True)

                processed[obs]['rebinned_extended'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                    input_data[obs]['step'],
                                    preparation[obs]['ratio'],
                                    blaze,
                                    processed['common_extended']['wave'],
                                    processed['common_extended']['step'],
                                    rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'],
                                    preserve_flux=preserve_flux)



            if norm_pams['normalize_transmission'] and norm_pams['normalization_model'] == 'polynomial':

                """ Continuum normalization preparatory steps:
                    1) exclusion of regions with planetary lines
                    2) exclusion of regions with stellar lines
                    3) Polynomial fit of selected regions
                        Boolean array initialized to all True values, fit is
                        performed on the extended region and then applied to the fit subset
                """
                processed['common_extended']['line_exclusion'] = (
                    processed['common_extended']['wave'] > 0.)

                """ Continuum normalization:
                    1) exclusion of regions with planetary lines, taking into
                        account the planetary RV semi-amplitude
                """
                for line_key, line_val in lines_dict['lines'].items():
                    line_extension = 1.2 * \
                        planet_dict['RV_semiamplitude'][0] * \
                            line_val / speed_of_light_km
                    processed['common_extended']['line_exclusion'] = processed['common_extended']['line_exclusion'] & (
                        np.abs(processed['common_extended']['wave']-line_val) > line_extension)

                """ Continuum normalization:
                    2) exclusion of regions with planetary lines, taking into
                        account the planetary RV semi-amplitude
                """

                try:
                    stellar_spectrum_rebinned = rebin_1d_to_1d(clv_rm_models['common']['wave'],
                                        clv_rm_models['common']['step'],
                                        clv_rm_models['common']['norm_convolved'],
                                        processed['common_extended']['wave'],
                                        processed['common_extended']['step'],
                                        preserve_flux=False)

                    stellar_spectrum_derivative = first_derivative(
                        processed['common_extended']['wave'], stellar_spectrum_rebinned)

                    cont_10perc = np.percentile(np.abs(stellar_spectrum_derivative), norm_pams['percentile_selection'])

                    processed['common_extended']['line_exclusion'] = processed['common_extended']['line_exclusion'] \
                            & (np.abs(stellar_spectrum_derivative) < cont_10perc) \
                            & (stellar_spectrum_rebinned > norm_pams['lower_threshold'])

                except KeyError:
                    print(
                        "No stellar synthetic spectrum from CLV models, some stellar lines may be included transmission normalization  ")

                for obs in lists['observations']:

                    selection = processed['common_extended']['line_exclusion'] & (
                        processed[obs]['rebinned_extended'] > np.std(processed[obs]['rebinned_extended']))

                    processed[obs]['norm_coeff'] = \
                        np.polynomial.chebyshev.chebfit(processed['common_extended']['wave'][selection],
                                                        processed[obs]['rebinned_extended'][selection],
                                                        norm_pams['spectra_poly_degree'])
                    processed[obs]['continuum'] = np.polynomial.chebyshev.chebval(
                        processed['common']['wave'], processed[obs]['norm_coeff'])
                    processed[obs]['normalized'] = processed[obs]['rebinned'] / \
                        processed[obs]['continuum']
                    processed[obs]['normalized_err'] = processed[obs]['rebinned_err'] / \
                        processed[obs]['continuum']

            elif norm_pams['normalize_transmission']  and (
                    norm_pams['normalization_model'] == 'savgol'
                    or norm_pams['normalization_model'] == 'savitzky-golay'):

                print('  Normalization using Savitzky-Golay filter')

                for obs in lists['observations']:

                    if norm_pams['normalize_rebinned']:

                        processed[obs]['continuum'] = savgol_filter(preparation[obs]['rebinned'],
                                                                window_length=norm_pams['window_length'],
                                                                polyorder=norm_pams['polyorder'],
                                                                mode=norm_pams['mode'],
                                                                cval=norm_pams['cval'])


                    else:
                        normalization_model = preparation[obs]['ratio'] * 0.00
                        for order in range(0, observational_pams['n_orders']):

                            normalization_model[order,:] = savgol_filter(preparation[obs]['ratio'][order,:],
                                                                window_length=norm_pams['window_length'],
                                                                polyorder=norm_pams['polyorder'],
                                                                mode=norm_pams['mode'],
                                                                cval=norm_pams['cval'])

                        processed[obs]['continuum'] = \
                            rebin_2d_to_1d(input_data[obs]['wave'],
                                            input_data[obs]['step'],
                                            normalization_model,
                                            blaze,
                                            processed['common']['wave'],
                                            processed['common']['step'],
                                            rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'],
                                            preserve_flux=preserve_flux)


                    processed[obs]['continuum_coeff'] = None
                    processed[obs]['normalized'] = processed[obs]['rebinned'] / \
                        processed[obs]['continuum']
                    processed[obs]['normalized_err'] = processed[obs]['rebinned_err'] / \
                        processed[obs]['continuum']


            else:
                for obs in lists['observations']:

                    processed[obs]['continuum_coeff'] = None
                    processed[obs]['continuum'] = np.ones_like(processed['common']['wave'])
                    processed[obs]['normalized'] = processed[obs]['rebinned'].copy()
                    processed[obs]['normalized_err'] = processed[obs]['rebinned_err'].copy()


            processed['common']['n_obs'] = len(lists['transit_full'])
            processed['common']['n_radius_grid'] = clv_rm_models['common']['n_radius_grid']
            processed['common']['radius_grid'] = clv_rm_models['common']['radius_grid']
            clv_rm_radius = clv_rm_models['common']['radius_grid']

            """ We are moving the values of interest from dictionaries to arrays
                in order to  speed up the MCMC
                1) clv_rm_grid: array with all the CLV models, as a function of the
                    radius of the planet
                2) time_from_transit: BJD_TDB - T0
                3) planet_RVsinusoid: Fractional RV of the planet (K=1) - from a meshgrid
            """
            clv_rm_grid = np.ones([processed['common']['n_radius_grid'],
                                    processed['common']['n_obs'],
                                    processed['common']['size']],
                                   dtype=np.double)

            time_from_transit = np.empty(
                processed['common']['n_obs'], dtype=np.double)

            transmission_spec = np.empty([processed['common']['n_obs'],
                                          processed['common']['size']],
                                         dtype=np.double)
            transmission_spec_err = np.empty([processed['common']['n_obs'],
                                              processed['common']['size']],
                                             dtype=np.double)

            for i_obs, obs in enumerate(lists['transit_full']):

                time_from_transit[i_obs] = observational_pams[obs]['BJD'] - \
                    observational_pams['time_of_transit']
                # planet_RVsinusoid[i_obs] = np.sin(2*np.pi /  planet_dict['period'][0] * time_from_transit[i_obs])

                transmission_spec[i_obs, :] = processed[obs]['normalized']
                transmission_spec_err[i_obs,
                                      :] = processed[obs]['normalized_err']

                if clv_rm_correction is False:
                    continue

                for i_r in range(0, processed['common']['n_radius_grid']):

                    """ CLV Synthetic models are in the Stellar Reference system,
                        so no shift is required """

                    clv_rm_grid[i_r, i_obs, :] = \
                        rebin_1d_to_1d(clv_rm_models['common']['wave'],
                                       clv_rm_models['common']['step'],
                                       clv_rm_models[obs]['clv_rm_model_convolved_normalized'][i_r, :],
                                       processed['common']['wave'],
                                       processed['common']['step'],
                                       preserve_flux=False)



                    # preserve_flux should be True or False?
                    # False if the spectra are already normalized


            #colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(
            #    lists, observational_pams)

            #fig = plt.figure(figsize=(12, 6))
            #gs = GridSpec(2, 2, width_ratios=[50, 1])

            #ax1 = plt.subplot(gs[0, 0])
            #ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
            #cbax1 = plt.subplot(gs[:, 1])

            #i_r = 0
            #for i_obs, obs in enumerate(lists['transit_full']):

            #    ax1.plot(processed['common']['wave'],
            #            clv_rm_grid[i_r, i_obs, :],
            #            color=colors_plot['mBJD'][obs], alpha=0.2)

            #i_r = processed['common']['n_radius_grid']-1
            #for i_obs, obs in enumerate(lists['transit_full']):

            #    ax2.plot(processed['common']['wave'],
            #            clv_rm_grid[i_r, i_obs, :],
            #            color=colors_plot['mBJD'][obs], alpha=0.2)

            #ax1.set_title(
            #    'Night: {0:s} \n CLV+RM correction, convolved and normalized '.format(night))
            #ax2.set_title('Out of transit')
            #ax2.set_xlabel('$\lambda$ [$\AA$]')

            #sm = plt.cm.ScalarMappable(
            #    cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
            #sm.set_array([])  # You have to set a dummy-array for this to work...
            #cbar = plt.colorbar(sm, cax=cbax1)
            #cbar.set_label('BJD - 2450000.0')
            #fig.subplots_adjust(wspace=0.05, hspace=0.4)
            #plt.show()

            #quit()


            remove_outliers = (np.abs(transmission_spec - 1.) > 0.5)
            transmission_spec[remove_outliers] = 1.0
            transmission_spec_err[remove_outliers] = 1.0
            wave_meshgrid, time_meshgrid = np.meshgrid(
                processed['common']['wave'], time_from_transit)
            planet_RVsinusoid = np.sin(
                2*np.pi / planet_dict['period'][0] * time_meshgrid)

            if jitter_flag:
                jitter_index = []
                n_jitter = 1
            else:
                jitter_index = None
                n_jitter = 0

            mcmc_data = {
                'observations': lists['transit_full'],
                'common_wave': processed['common']['wave'],
                'common_step': processed['common']['step'],
                'clv_rm_grid': clv_rm_grid,
                'transmission_spec': transmission_spec,
                'transmission_spec_err': transmission_spec_err,
                'wave_meshgrid': wave_meshgrid,
                'time_meshgrid': time_meshgrid,
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
        print("transmission_binned_mcmc           ")

        try:
            results_dict = load_from_cpickle(subroutine_name+'_'+sampler_name+'_results', config_in['output'], night, lines_label, it_string)
            print("   Transmission MCMC analysis for lines {0:s}, night: {1:s}  already performed".format(
                lines_label, night))

            pams_dict = results_dict['pams_dict']
            chain_med = results_dict['chain_med']
            boundaries = results_dict['boundaries']
            start_average = np.average(results_dict['point_start'], axis=0)
            ndim = results_dict['ndim']
            med_lines_model = results_dict['results']['lines_model']

            if 'derived' in results_dict:
                recompute_derived = False
            else:
                recompute_derived = True
                results_dict['derived'] = {}

            # TODO improve output
            print('   *** sampler output ')

            for key, val in pams_dict.items():
                print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}]) (start: {7:9f})'.format(key, val,
                                                                                                        chain_med[val,0],
                                                                                                        chain_med[val,2],
                                                                                                        chain_med[val,1],
                                                                                                        boundaries[val, 0],
                                                                                                        boundaries[val, 1],
                                                                                                        start_average[val])
                                                                                                        )
                if recompute_derived and key[-8:]=='contrast':
                    key_name = key[:-8] + 'Rh'
                    sample_size = len(results_dict['flat_chain'][:,val])
                    planet_ratio_sample = np.random.normal(planet_dict['radius_ratio'][0],planet_dict['radius_ratio'][1],size=sample_size)
                    results_dict['derived'][key_name] = {}
                    results_dict['derived'][key_name]['flat_chain'] = np.sqrt(results_dict['flat_chain'][:,val]/planet_ratio_sample**2 + 1.)
                    results_dict['derived'][key_name]['chain_med'] = compute_value_sigma(results_dict['derived'][key_name]['flat_chain'])

                #    R(h) = np.sqrt(1+h/delta)

                #    print(key[-8:], key[:3])

            print('   *** derived output ')
            for key, val in results_dict['derived'].items():
                chain_med = results_dict['derived'][key]['chain_med']
                print('{0:24s}         {1:12f}   {2:12f}  {3:12f} (15-84 p)'.format(key,
                                                                                     chain_med[0],
                                                                                     chain_med[2],
                                                                                     chain_med[1]))



            continue



        except FileNotFoundError:
            print()

        # getting fit parameters
        lines_center, pams_dict, pams_list, boundaries, theta_start = define_theta_array(
            model_case, lines_dict, planet_dict, n_jitter, allow_emission=allow_emission)
        ndim = len(theta_start)

        if pyde_flag:
            ngen = sampler_pams.get('n_gen', 64000)
        else:
            ngen = 0
        nwalkers_mult = sampler_pams.get('n_walkers_mult', 2)
        nwalkers = sampler_pams.get('n_walkers', nwalkers_mult * ndim)
        nthin = sampler_pams.get('n_thin', 50)
        nsteps = sampler_pams.get('n_steps', 20000)
        nburnin = sampler_pams.get('n_burnin', 10000)
        ncpus = sampler_pams.get('n_cpus', 4)
        ndata = np.size(wave_meshgrid)

        if pams_dict.get('rp_factor', False):
            pam_id = pams_dict['rp_factor']
            boundaries[pam_id, :] = [clv_rm_radius[0], clv_rm_radius[-1]]


        print()
        print('      PyDE + emcee parameters')
        print('         n_dim:                        {0:9.0f}'.format(ndim))
        if pyde_flag:
            print('         n_gen:      (default: 64000)  {0:9.0f}'.format(ngen))
        else:
            print('         no PyDE optimization, MCMC will start from default values')
        print(
            '         n_walkers:  (default: 2*ndim) {0:9.0f}'.format(nwalkers))
        print('         n_steps:    (default: 20000)  {0:9.0f}'.format(nsteps))
        print(
            '         n_burnin:   (default: 10000)  {0:9.0f}'.format(nburnin))
        print('         n_thin:     (default: 50)     {0:9.0f}'.format(nthin))

        population, sampler_chain, sampler_lnprobability, point_start = emcee_lines_fit_functions(
            model_case,
            wave_meshgrid,
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
                         wave_meshgrid,
                         clv_rm_radius,
                         clv_rm_grid,
                         planet_RVsinusoid,
                         lines_center,
                         jitter_index)

        map_lines_model, map_clv_model, map_lines_array, map_planet_K, map_planet_R, map_jitter = \
            return_model(model_case,
                         chain_MAP,
                         wave_meshgrid,
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

            # ================================================================
            # RV shift from the Observer RF to the Planet RF
            # Supports both circular and eccentric orbits using compute_planet_RV
            #STRONG ASSUMPTIONS:
            # - there is only the transiting planet in the system
            # - linear approximation or the orbit near the transit event 
            # ================================================================
            # Extract orbital parameters (set during prepare_datasets)
            orb_pams = observational_pams.get('orbital_parameters', {
                'eccentricity': 0.0,
                'omega_rad': np.pi / 2.0,
                'reference_time': observational_pams['time_of_transit']
            })
            
            # Compute planet RV using MCMC-derived K value (results)
            rv_planet_results = compute_planet_RV(
                BJD=observational_pams[obs]['BJD'],
                time_of_transit=observational_pams['time_of_transit'],
                period=planet_dict['period'][0],
                K_planet=results_dict['results']['planet_K'],
                e0=orb_pams['eccentricity'],
                omega0=orb_pams['omega_rad'],
                reference_time=orb_pams['reference_time'],
                output_unit='km/s'
            )
            
            # Compute planet RV using MCMC-derived K value (results_MAP)
            rv_planet_map = compute_planet_RV(
                BJD=observational_pams[obs]['BJD'],
                time_of_transit=observational_pams['time_of_transit'],
                period=planet_dict['period'][0],
                K_planet=results_dict['results_MAP']['planet_K'],
                e0=orb_pams['eccentricity'],
                omega0=orb_pams['omega_rad'],
                reference_time=orb_pams['reference_time'],
                output_unit='km/s'
            )
            
            # RV shift from Observer RF to Planet RF
            results_dict['results']['observational_pams'][obs]['rv_shift_ORF2PRF'] = \
                observational_pams[obs]['BERV'] \
                - observational_pams['RV_star']['RV_systemic'] \
                - rv_planet_results

            results_dict['results_MAP']['observational_pams'][obs]['rv_shift_ORF2PRF'] = \
                observational_pams[obs]['BERV'] \
                - observational_pams['RV_star']['RV_systemic'] \
                - rv_planet_map

            # RV shift from Stellar Rest Frame to Planetary Rest Frame
            results_dict['results']['observational_pams'][obs]['rv_shift_SRF2PRF'] = \
                + observational_pams[obs]['RV_bjdshift'] \
                - rv_planet_results

            results_dict['results_MAP']['observational_pams'][obs]['rv_shift_SRF2PRF'] = \
                + observational_pams[obs]['RV_bjdshift'] \
                - rv_planet_map


        results_dict['derived'] = {}

        # TODO improve output
        print('   *** sampler output ')

        start_average = np.average(results_dict['point_start'], axis=0)
        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}]) (start: {7:9f})'.format(key, val,
                                                                                                    chain_med[val,0],
                                                                                                    chain_med[val,2],
                                                                                                    chain_med[val,1],
                                                                                                    boundaries[val, 0],
                                                                                                    boundaries[val, 1],
                                                                                                    start_average[val])
                                                                                                    )


            if key[-8:]=='contrast':
                key_name = key[:-8] + 'Rh'
                sample_size = len(results_dict['flat_chain'][:,val])
                planet_ratio_sample = np.random.normal(planet_dict['radius_ratio'][0],planet_dict['radius_ratio'][1],size=sample_size)
                results_dict['derived'][key_name] = {}
                results_dict['derived'][key_name]['flat_chain'] = np.sqrt(results_dict['flat_chain'][:,val]/planet_ratio_sample**2 + 1.)
                results_dict['derived'][key_name]['chain_med'] = compute_value_sigma(results_dict['derived'][key_name]['flat_chain'])


        print('   *** derived output ')
        for key, val in results_dict['derived'].items():
            chain_med = results_dict['derived'][key]['chain_med']
            print('{0:24s}         {1:12f}   {2:12f}  {3:12f} (15-84 p)'.format(key,
                                                                                    chain_med[0],
                                                                                    chain_med[2],
                                                                                    chain_med[1]))


        save_to_cpickle(subroutine_name+'_'+sampler_name+'_results',
                        results_dict, config_in['output'], night, lines_label, it_string)

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
        all_mcmc_data = load_from_cpickle(subroutine_name+'_data', config_in['output'], night='', lines=lines_label, it_string=it_string)

        all_clv_rm_radius = all_mcmc_data['clv_rm_radius']
        all_clv_rm_grid = all_mcmc_data['clv_rm_grid']
        all_transmission_spec = all_mcmc_data['transmission_spec']
        all_transmission_spec_err = all_mcmc_data['transmission_spec_err']
        all_wave_meshgrid = all_mcmc_data['wave_meshgrid']
        all_time_meshgrid = all_mcmc_data['time_meshgrid']
        all_planet_RVsinusoid = all_mcmc_data['planet_RVsinusoid']
        all_observations = all_mcmc_data['observations']
        all_n_obs = all_mcmc_data['n_obs']
        all_n_radius_grid = all_mcmc_data['n_radius_grid']
        all_jitter_index = all_mcmc_data['jitter_index']
        n_jitter = all_mcmc_data['n_jitter']

    except:

        n_jitter = 0

        for night in night_dict:

            mcmc_data = load_from_cpickle(subroutine_name+'_data', config_in['output'], night, lines_label, it_string=it_string)

            try:
                # Building the arrays for the full analysis

                all_clv_rm_grid = np.concatenate(
                    (all_clv_rm_grid, mcmc_data['clv_rm_grid']), axis=1)
                all_transmission_spec = np.concatenate(
                    (all_transmission_spec, mcmc_data['transmission_spec']))
                all_transmission_spec_err = np.concatenate(
                    (all_transmission_spec_err, mcmc_data['transmission_spec_err']))
                all_wave_meshgrid = np.concatenate(
                    (all_wave_meshgrid, mcmc_data['wave_meshgrid']))
                all_time_meshgrid = np.concatenate(
                    (all_time_meshgrid, mcmc_data['time_meshgrid']))
                all_planet_RVsinusoid = np.concatenate(
                    (all_planet_RVsinusoid, mcmc_data['planet_RVsinusoid']))
                all_observations = np.concatenate(
                    (all_observations, mcmc_data['observations']))
                all_n_obs += mcmc_data['n_obs']

                if jitter_flag:
                    all_jitter_index = np.concatenate(
                        (all_jitter_index, n_jitter*np.ones(np.shape(mcmc_data['wave_meshgrid']), dtype=np.int16)))
                    n_jitter += 1

            except NameError:
                """ This error is expected when retrieving the data of the first night"""
                all_clv_rm_radius = mcmc_data['clv_rm_radius']
                all_clv_rm_grid = mcmc_data['clv_rm_grid']
                all_transmission_spec = mcmc_data['transmission_spec']
                all_transmission_spec_err = mcmc_data['transmission_spec_err']
                all_wave_meshgrid = mcmc_data['wave_meshgrid']
                all_time_meshgrid = mcmc_data['time_meshgrid']
                all_planet_RVsinusoid = mcmc_data['planet_RVsinusoid']
                all_observations = mcmc_data['observations']
                all_n_obs = mcmc_data['n_obs']
                all_n_radius_grid = mcmc_data['n_radius_grid']

                if jitter_flag:
                    all_jitter_index = n_jitter * \
                        np.ones(
                            np.shape(mcmc_data['wave_meshgrid']), dtype=np.int16)
                    n_jitter += 1
                else:
                    all_jitter_index = None

        all_mcmc_data = {
            'observations': all_observations,
            'clv_rm_grid': all_clv_rm_grid,
            'transmission_spec': all_transmission_spec,
            'transmission_spec_err': all_transmission_spec_err,
            'wave_meshgrid': all_wave_meshgrid,
            'time_meshgrid': all_time_meshgrid,
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
        results_dict = load_from_cpickle(subroutine_name+ '_'+ sampler_name+'_results',
                                         config_in['output'], night='', lines=lines_label, it_string=it_string)
        print("   Transmission MCMC analysis for lines {0:s} already performed ".format(
                lines_label))

        pams_dict = results_dict['pams_dict']
        chain_med = results_dict['chain_med']
        boundaries = results_dict['boundaries']
        ndim = results_dict['ndim']
        start_average = np.average(results_dict['point_start'], axis=0)

        if 'derived' in results_dict:
            recompute_derived = False
        else:
            recompute_derived = True
            results_dict['derived'] = {}

        # TODO improve output
        print('   *** sampler output ')

        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}]) (start: {7:9f})'.format(key, val,
                                                                                                    chain_med[val,0],
                                                                                                    chain_med[val,2],
                                                                                                    chain_med[val,1],
                                                                                                    boundaries[val, 0],
                                                                                                    boundaries[val, 1],
                                                                                                    start_average[val])
                                                                                                    )

            if recompute_derived and key[-8:]=='contrast':
                key_name = key[:-8] + 'Rh'
                sample_size = len(results_dict['flat_chain'][:,val])
                planet_ratio_sample = np.random.normal(planet_dict['radius_ratio'][0],planet_dict['radius_ratio'][1],size=sample_size)
                results_dict['derived'][key_name] = {}
                results_dict['derived'][key_name]['flat_chain'] = np.sqrt(results_dict['flat_chain'][:,val]/planet_ratio_sample**2 + 1.)
                results_dict['derived'][key_name]['chain_med'] = compute_value_sigma(results_dict['derived'][key_name]['flat_chain'])
            #    R(h) = np.sqrt(1+h/delta)

                #    print(key[-8:], key[:3])

        print('   *** derived output ')
        for key, val in results_dict['derived'].items():
            chain_med = results_dict['derived'][key]['chain_med']
            print('{0:24s}         {1:12f}   {2:12f}  {3:12f} (15-84 p)'.format(key,
                                                                                    chain_med[0],
                                                                                    chain_med[2],
                                                                                    chain_med[1]))





    except FileNotFoundError:

        lines_center, pams_dict, pams_list, boundaries, theta_start = define_theta_array(
            model_case, lines_dict, planet_dict, n_jitter, allow_emission=allow_emission)
        ndim = len(theta_start)
        ngen = sampler_pams.get('n_gen', 64000)
        nwalkers_mult = sampler_pams.get('n_walkers_mult', 2)
        nwalkers = sampler_pams.get('n_walkers', nwalkers_mult * ndim)
        nthin = sampler_pams.get('n_thin', 50)
        nsteps = sampler_pams.get('n_steps', 20000)
        nburnin = sampler_pams.get('n_burnin', 10000)
        ncpus = sampler_pams.get('n_cpus', 4)
        ndata = np.size(all_wave_meshgrid)

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
            all_wave_meshgrid,
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
                         all_wave_meshgrid,
                         all_clv_rm_radius,
                         all_clv_rm_grid,
                         all_planet_RVsinusoid,
                         lines_center,
                         all_jitter_index)

        map_lines_model, map_clv_model, map_lines_array, map_planet_K, map_planet_R, map_jitter = \
            return_model(model_case,
                         chain_MAP,
                         all_wave_meshgrid,
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
            'theta_start': theta_start
            #'BIC': BIC,
            #'BIC_map': BIC_map
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



        start_average = np.average(results_dict['point_start'], axis=0)

        results_dict['derived'] = {}

        # TODO improve output
        print('   *** sampler output ')
        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}]) (start: {7:9f})'.format(key, val,
                                                                                                    chain_med[val,0],
                                                                                                    chain_med[val,2],
                                                                                                    chain_med[val,1],
                                                                                                    boundaries[val, 0],
                                                                                                    boundaries[val, 1],
                                                                                                    start_average[val])
                                                                                                    )
            if key[-8:]=='contrast':
                key_name = key[:-8] + 'Rh'
                sample_size = len(results_dict['flat_chain'][:,val])
                print(sample_size)
                planet_ratio_sample = np.random.normal(planet_dict['radius_ratio'][0],planet_dict['radius_ratio'][1],size=sample_size)
                results_dict['derived'][key_name] = {}
                results_dict['derived'][key_name]['flat_chain'] = np.sqrt(results_dict['flat_chain'][:,val]/planet_ratio_sample**2 + 1.)
                results_dict['derived'][key_name]['chain_med'] = compute_value_sigma(results_dict['derived'][key_name]['flat_chain'])


        print('   *** derived output ')
        for key, val in results_dict['derived'].items():
            chain_med = results_dict['derived'][key]['chain_med']
            print('{0:24s}         {1:12f}   {2:12f}  {3:12f} (15-84 p)'.format(key,
                                                                                    chain_med[0],
                                                                                    chain_med[2],
                                                                                    chain_med[1]))

        save_to_cpickle(subroutine_name +'_'+sampler_name+'_results',
                        results_dict, config_in['output'], night='', lines=lines_label, it_string=it_string)


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



def plot_transmission_binned_mcmc(config_in, lines_label, night_input='', reference='planetRF', pca_iteration=-1):

    night_dict = from_config_get_nights(config_in)
    planet_dict = from_config_get_planet(config_in)
    star_dict = from_config_get_star(config_in)
    clv_rm_dict = from_config_get_clv_rm(config_in)

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]
    sampler_pams = lines_dict['sampler_parameters']
    sampler_name = sampler_pams.get('sampler_name', 'emcee')

    if night_input == '':
        night_list = ['']
    else:
        night_list = np.atleast_1d(night_input)

    os.system('mkdir -p plots')


    # Workaround to check if the transmission spectrum has been obtained through PCA iterations
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


    for night in night_list:

        results_dict = load_from_cpickle(subroutine_name+'_'+sampler_name+'_results', config_in['output'], night, lines_label, it_string)
        print("   Transmission MCMC analysis for lines {0:s}, night: {1:s}  already performed".format(
            lines_label, night))

        if night == '':
            chains_dir = 'plots/mcmc_binned_chains_full/'
        else:
            chains_dir = 'plots/mcmc_binned_chains_' + night + '/'
        os.system('mkdir -p ' + chains_dir)


        pams_dict = results_dict['pams_dict']
        chain_med = results_dict['chain_med']
        lnprob_med = results_dict['lnprob_med']
        boundaries = results_dict['boundaries']
        flat_chain = results_dict['flat_chain']
        flat_lnprob = results_dict['flat_lnprob']
        nthin = results_dict['nthin']
        nsteps = results_dict['nsteps']
        nburnin = results_dict['nburnin']

        sampler_chain = results_dict['sampler_chain']

        start_average = np.average(results_dict['point_start'], axis=0)
        ndim = results_dict['ndim']
        med_lines_model = results_dict['results']['lines_model']

        if 'derived' in results_dict:
            recompute_derived = False
        else:
            recompute_derived = True
            results_dict['derived'] = {}

        # TODO improve output
        print('   *** sampler output (plotting the chains)')

        sample_size = np.size(flat_chain, axis=0)
        dimen_size = np.size(flat_chain, axis=1)
        corner_plot = {
            'samples': np.zeros([sample_size, dimen_size + 1]),
            'labels': [],
            'truths': [],
            'start': [],
        }

        i_corner = 0


        for key, val in pams_dict.items():
            print('{0:24s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}]) (start: {7:9f})'.format(key, val,
                                                                                                    chain_med[val,0],
                                                                                                    chain_med[val,2],
                                                                                                    chain_med[val,1],
                                                                                                    boundaries[val, 0],
                                                                                                    boundaries[val, 1],
                                                                                                    start_average[val])
                                                                                                    )
            if recompute_derived and key[-8:]=='contrast':
                key_name = key[:-8] + 'Rh'
                planet_ratio_sample = np.random.normal(planet_dict['radius_ratio'][0],planet_dict['radius_ratio'][1],size=sample_size)
                results_dict['derived'][key_name] = {}
                results_dict['derived'][key_name]['flat_chain'] = np.sqrt(results_dict['flat_chain'][:,val]/planet_ratio_sample**2 + 1.)
                results_dict['derived'][key_name]['chain_med'] = compute_value_sigma(results_dict['derived'][key_name]['flat_chain'])


            corner_plot['samples'][:, i_corner] = flat_chain[:, val]
            corner_plot['labels'].append(re.sub('_', ' ', key))
            corner_plot['truths'].append(chain_med[val, 0])
            corner_plot['start'].append(start_average[val])
            i_corner += 1


            file_name = chains_dir + repr(val) + '.png'
            fig = plt.figure(figsize=(12, 12))
            plt.title(key)
            plt.plot(sampler_chain[:, :, val].T, '-', alpha=0.5)
            plt.axvline(nburnin / nthin, c='r')
            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)



        corner_plot['samples'][:, -1] = flat_lnprob[:]
        corner_plot['labels'].append('ln-prob')
        corner_plot['truths'].append(lnprob_med[0])
        corner_plot['start'].append(None)

            #    R(h) = np.sqrt(1+h/delta)

            #    print(key[-8:], key[:3])

        print('   *** derived output ')
        for key, val in results_dict['derived'].items():
            chain_med = results_dict['derived'][key]['chain_med']
            print('{0:24s}         {1:12f}   {2:12f}  {3:12f} (15-84 p)'.format(key,
                                                                                    chain_med[0],
                                                                                    chain_med[2],
                                                                                    chain_med[1]))
        print('   *** corner plot using pyGTC output ')

        filename_rad = subroutine_name + '_' + reference + '_cornerplot'
        output_file = get_filename(filename_rad, config_in['output'], night=night, lines=lines_label, it_string=it_string, extension='.pdf')
        print('   *** filename: ', output_file)

        try:
            GTC = pygtc.plotGTC(chains=corner_plot['samples'],
                                paramNames=corner_plot['labels'],
                                truths=[corner_plot['truths'],corner_plot['start']],
                                GaussianConfLevels=True,
                                nConfidenceLevels=3,
                                figureSize=12,
                                labelRotation= (True,True),
                                plotName='plots/'+output_file)
        except:
            GTC = pygtc.plotGTC(chains=corner_plot['samples'],
                                paramNames=corner_plot['labels'],
                                truths=[corner_plot['truths'],corner_plot['start']],
                                GaussianConfLevels=True,
                                nConfidenceLevels=2,
                                figureSize=12,
                                labelRotation= (True,True),
                                plotName='plots/'+output_file)


        GTC = None


        continue


def plot_transmission_binned_mcmc_deprecated(config_in, lines_label, night_input=''):
    night_dict = from_config_get_nights(config_in)

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the analysis"""
        try:
            clv_rm_corrected = load_from_cpickle(subroutine_name, config_in['output'], night)

            mcmc_data = load_from_cpickle(subroutine_name + '_data', config_in['output'], night, lines_label)

            clv_rm_radius = mcmc_data['clv_rm_radius']
            clv_rm_grid = mcmc_data['clv_rm_grid']
            transmission_spec = mcmc_data['transmission_spec']
            transmission_spec_err = mcmc_data['transmission_spec_err']
            wave_meshgrid = mcmc_data['wave_meshgrid']
            time_meshgrid = mcmc_data['time_meshgrid']
            planet_RVsinusoid = mcmc_data['planet_RVsinusoid']

        except:
            print("No transmission spectrum results, no plots")
            print()
            continue

        """ Creation of the color array, based on the BJD of the observations
        """
        bjd = []
        am = []

        for obs in lists['observations']:
            bjd.append(clv_rm_corrected[obs]['BJD'] - 2450000.0)
            am.append(clv_rm_corrected[obs]['AIRMASS'])



        from matplotlib.colors import BoundaryNorm
        from matplotlib.ticker import MaxNLocator
        cmap = plt.get_cmap('coolwarm')

        plot_data = transmission_spec.copy()


        from SLOPpy.subroutines.math_functions import interpolate2d_grid_nocheck
        plot_data = interpolate2d_grid_nocheck(1.000, clv_rm_radius, clv_rm_grid)
        #clv_model = interpolate2d_grid_nocheck(1.000, clv_rm_radius, clv_rm_grid)
        vmin = plot_data.min()
        vmax =  plot_data.max()

        levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plt.figure(figsize=(15, 10))

        PCF = plt.contourf(wave_meshgrid, time_meshgrid, plot_data, levels=levels, cmap=cmap)
        cbar = plt.colorbar(PCF)
        cbar.ax.set_ylabel('Intensity')
        plt.show()

        levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plt.figure(figsize=(15, 10))

        PCM = plt.pcolormesh(wave_meshgrid, time_meshgrid, plot_data,
                             vmin=vmin, vmax=vmax, cmap=cmap)
        cbar = plt.colorbar(PCM)
        cbar.ax.set_ylabel('Intensity')
        plt.show()




        plot_data = transmission_spec.copy()


        plot_data = transmission_spec.copy()
        #clv_model = interpolate2d_grid_nocheck(1.000, clv_rm_radius, clv_rm_grid)
        vmin = plot_data.min()
        vmax =  plot_data.max()

        levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plt.figure(figsize=(15, 10))

        PCF = plt.contourf(wave_meshgrid, time_meshgrid, plot_data, levels=levels, cmap=cmap)
        cbar = plt.colorbar(PCF)
        cbar.ax.set_ylabel('Intensity')
        plt.show()

        levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plt.figure(figsize=(15, 10))

        PCM = plt.pcolormesh(wave_meshgrid, time_meshgrid, plot_data,
                             vmin=vmin, vmax=vmax, cmap=cmap)
        cbar = plt.colorbar(PCM)
        cbar.ax.set_ylabel('Intensity')
        plt.show()
