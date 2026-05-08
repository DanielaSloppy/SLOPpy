from __future__ import print_function, division
from SLOPpy.instruments.get_data import *
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.kepler_exo import compute_planet_RV

__all__ = ["prepare_datasets", "plot_dataset"]


def _check_wavelength_rescaling(wave_rescaling, wave_observations):
    if wave_rescaling[0] < wave_observations[0] or \
            wave_rescaling[1] > wave_observations[1]:

        warnings.warn("Valid wavelength rescaling window must be between {0:8.2f} and {1:8.2f}".format(
            wave_observations[0], wave_observations[1]))

        return False
    else:
        return True


def _check_coadd_in_shared_data(shared_data, wavelength_range):
    if 'coadd' not in shared_data:
        shared_data['coadd'] = {
            'wavelength_range': wavelength_range[:]
        }
        shared_data['binned'] = {
            'wavelength_range': wavelength_range[:]
        }
    else:
        shared_data['coadd']['wavelength_range'][0] = min(shared_data['coadd']['wavelength_range'][0],
                                                          wavelength_range[0])

        shared_data['coadd']['wavelength_range'][1] = max(shared_data['coadd']['wavelength_range'][1],
                                                          wavelength_range[1])

        shared_data['binned']['wavelength_range'][0] = min(shared_data['binned']['wavelength_range'][0],
                                                           wavelength_range[0])

        shared_data['binned']['wavelength_range'][1] = max(shared_data['binned']['wavelength_range'][1],
                                                           wavelength_range[1])

    return shared_data


def prepare_datasets(config_in):
    """
    Prepare observation datasets for the SLOPpy transmission spectroscopy pipeline.
    
    This function reads FITS files and associated data, processes them, and saves
    them as pickle objects for subsequent pipeline steps. It handles multiple
    observation nights and instruments, computing observational parameters including
    radial velocity shifts for different reference frames.
    
    The function supports both circular and eccentric planetary orbits. For eccentric
    orbits, the full Keplerian solution is used to compute planet radial velocities.
    
    Parameters
    ----------
    config_in : dict
        Configuration dictionary containing all pipeline parameters, including:
        - output: str - Output directory/prefix for pickle files
        - nights: dict - Dictionary of observation nights with their parameters
        - instruments: dict - Dictionary of instrument configurations
        - planet: dict - Planetary parameters (period, RV_semiamplitude, etc.)
        - star: dict - Stellar parameters
        - master-out: dict - Master-out spectrum parameters
        
        For eccentric orbits, planet dict should include:
        - orbit: 'eccentric' (default is 'circular')
        - eccentricity: [value, error]
        - omega: [value, error] in degrees
    
    Returns
    -------
    None
        Results are saved to pickle files:
        - lists: observation file lists
        - input_dataset_fibA: fiber A spectral data
        - input_dataset_s1d_fibA: fiber A 1D spectra
        - calibration_fibA: calibration data for fiber A
        - observational_pams: computed observational parameters including RV shifts
        - shared: shared wavelength grid and coadd parameters
        
        Optional files (if fiber B data available):
        - input_dataset_fibB: fiber B spectral data
        - calibration_fibB: calibration data for fiber B
    
    Notes
    -----
    This function performs the following main tasks:
    
    1. Reads observation file lists for each night
    2. Loads spectral data from FITS files
    3. Handles negative flux values and estimates noise floor
    4. Determines transit in/out observations based on transit duration
    5. Computes RV shifts between reference frames:
       - Observer Reference Frame (ORF)
       - Solar Barycentric Reference Frame (BRF)
       - Stellar Reference Frame (SRF)
       - Planetary Reference Frame (PRF)
    6. Sets up common wavelength grids for coadding spectra
    
    For existing datasets, the function loads previously computed parameters
    and only updates the transit lists if necessary.
    
    See Also
    --------
    _get_observational_parameters : Compute RV shifts for each observation
    compute_planet_RV : Compute planet RV for circular/eccentric orbits
    """

    """ config_dictionary: dictionary with all the configuration parameters from config_in
        lists_dictionary: for each night, the list of files
    """

    night_dict = from_config_get_nights(config_in)

    instrument_dict = from_config_get_instrument(config_in)

    #try:
    #    shared_data = load_from_cpickle('shared', config_in['output'])
    #    loaded_shared_data = True
    #except:
    #    shared_data = {}
    #    loaded_shared_data = False

    if check_existence_cpickle('shared', config_in['output']):
        loaded_shared_data = True
    else:
        shared_data = {}
        loaded_shared_data = False

    pass_wavelength_rescaling = True
    pass_wavelength_master_out = True

    for night in night_dict:

        print('Processing data for night: ', night)
        print()

        """ List files are supposed to be in the same directory of the yaml file,
            NOT on the archive directory: in this way it is possible to try different
            combinations of nights and files without making a mess in the archive """
        files_list, files_transit_out, files_transit_in, files_transit_full, files_telluric, files_star_telluric = get_filelists(
            night_dict[night])

        """ fix the potential problem of multiple observations taken during the same night"""
        night_name = night_dict[night].get('night', night)

        lists_dictionary = {
            'observations': files_list,
            'star_telluric': files_star_telluric,
            'n_observations': len(files_list),
        }

        try:
            lists_dictionary['transit_out'] = files_transit_out
            lists_dictionary['transit_in'] = files_transit_in
            lists_dictionary['transit_full'] = files_transit_full
            lists_dictionary['n_transit_out'] = len(files_transit_out)
            lists_dictionary['n_transit_in'] = len(files_transit_in)
            lists_dictionary['n_transit_full'] = len(files_transit_full)
            lists_dictionary['telluric'] = files_telluric
            lists_dictionary['n_tellurics'] = len(files_telluric)
            write_transit_list = False
        except:
            print('  Input lists for transit in/out not found, proceeding to automatic selection and writing ')
            print()
            write_transit_list = True

        """ Retrieval on instrument characteristics """
        instrument = night_dict[night]['instrument']
        mask = night_dict[night].get('mask', 'undefined')
        archive_dir = instrument_dict[instrument]['data_archive']
        order_selection = instrument_dict[instrument]['orders']
        wavelength_rescaling = instrument_dict[instrument]['wavelength_rescaling']
        time_of_transit = night_dict[night]['time_of_transit']

        planet_dict = from_config_get_planet(config_in)
        star_dict = from_config_get_star(config_in)

        print("  # observations: ", lists_dictionary['n_observations'])
        try:
            print("  # out-transit obs: ", lists_dictionary['n_transit_out'])
            print("  # in-transit obs: ", lists_dictionary['n_transit_in'])
        except:
            pass
        print("  Instrument: ", instrument)
        print("  Archive DIR: ", archive_dir)
        print("  Night: ", night)
        print("  Mask:  ", mask)
        print("  WL rescaling: ", wavelength_rescaling)
        print()

        try:

            if not check_existence_cpickle('input_dataset_fibA', config_in['output'], night):
                raise ValueError()

            #pass_wavelength_rescaling = _check_wavelength_rescaling(wavelength_rescaling,
            #                                                        observations_A['coadd']['wavelength_range']
            #                                                        )

            #shared_data = _check_coadd_in_shared_data(shared_data,
            #                                          observations_A['coadd']['wavelength_range'])
            print("  Input data for night {0:s} successfully retrieved".format(night))

            if write_transit_list:
                observations_A = load_from_cpickle('input_dataset_fibA', config_in['output'], night)

                lists_dictionary = _write_transit_list(observations_A,
                                                       lists_dictionary,
                                                       night_dict[night],
                                                       planet_dict)

                save_to_cpickle('lists', lists_dictionary, config_in['output'], night)
                print("  List rewritten, be careful however that you may incur in extra problems if they have changed")

            try:
                observational_parameters = load_from_cpickle('observational_pams', config_in['output'], night)
            except:
                observations_A = load_from_cpickle('input_dataset_fibA', config_in['output'], night)
                observational_parameters = _get_observational_parameters(observations_A,
                                                                         lists_dictionary,
                                                                         night_dict[night],
                                                                         instrument_dict[instrument],
                                                                         star_dict,
                                                                         planet_dict)
                save_to_cpickle('observational_pams', observational_parameters, config_in['output'], night)
                print("  New observational parameters loaded successfully")

            for key_name, key_val in observational_parameters['RV_star'].items():
                print("   RV star  {0:s}: {1}".format(key_name, key_val))
            print()

            # Ensure shared_data wavelength range is populated even when
            # input_dataset_fibA is already cached (shared.p may be missing).
            if not loaded_shared_data:
                _obs_tmp = load_from_cpickle('input_dataset_fibA', config_in['output'], night)
                shared_data = _check_coadd_in_shared_data(
                    shared_data, _obs_tmp['coadd']['wavelength_range'])
            continue

        except ValueError:
            pass

        observations_A = {}
        observations_s1d_A = {}
        observations_B = {}
        calib_data_A = {}
        calib_data_B = {}

        for obs in lists_dictionary['observations']:

            print("  Reading ", obs, " associated files")

            observations_A[obs], observations_s1d_A[obs] = \
                get_input_data(instrument, night_dict[night], archive_dir + night_name, obs, mask, skip_s1d=False,
                               order_selection=order_selection)   

            #""" Zero or negative values are identified, flagged and substituted with another value """
            #replacement = 0.01
            #observations_A[obs]['null'] = (observations_A[obs]['e2ds'] <= replacement)
            #observations_A[obs]['e2ds'][observations_A[obs]['null']] = replacement

            """ Negative values are just statistical noise around the null flux points,
            removing them would bias the flux level of the sky towards
            higher values

            We proceed in this way:
            -  identify the negative values and make a statistics of their
                average value
            - if in relevant number, we assume that the median of their absolute
                values corresponds to the noise floor
            - we add the noise floor to the error estimate 
            """

            observations_A[obs]['null'] = (observations_A[obs]['e2ds'] <= 0.0)
            if (np.sum(observations_A[obs]['null']) > 30):
                observations_A[obs]['noise_floor'] = np.median(np.abs(
                    observations_A[obs]['e2ds'][observations_A[obs]['null']]))
            else:
                if observations_A[obs].get('absolute_flux', True):
                    observations_A[obs]['noise_floor'] = 1.0000
                else:
                    observations_A[obs]['noise_floor'] = 0.00001


            observations_A[obs]['e2ds_err'] = np.sqrt(observations_A[obs]['e2ds_err']**2 + observations_A[obs]['noise_floor']**2)
            #observations_A[obs]['e2ds_err'] = np.sqrt(observations_A[obs]['e2ds'])

            if 'n_orders' not in observations_A or 'n_pixels' not in observations_A:
                observations_A['n_orders'] = observations_A[obs]['n_orders']
                observations_A['n_pixels'] = observations_A[obs]['n_pixels']
                calib_data_A = get_calib_data(instrument, night_dict[night], archive_dir + night_name, obs,
                                              order_selection=order_selection)

            """ Updating info on shared data """
            if 'coadd' not in observations_A:
                observations_A['coadd'] = {
                    'wavelength_range': [np.min(observations_A[obs]['wave'][0, :]),
                                         np.max(observations_A[obs]['wave'][-1, :])]
                }

            else:
                observations_A['coadd']['wavelength_range'][0] = min(observations_A['coadd']['wavelength_range'][0],
                                                                     np.min(observations_A[obs]['wave'][0, :]))

                observations_A['coadd']['wavelength_range'][1] = max(observations_A['coadd']['wavelength_range'][1],
                                                                     np.max(observations_A[obs]['wave'][-1, :]))

            """ Reading the fiber B counter part
                If the target has been observed in ThAr or FP mode, fiber B data will not be accessible
            """
            has_fiber_B = False
            try:

                observations_B[obs], _ = get_input_data(instrument, night_dict[night], archive_dir + night_name, obs, mask,
                                                        fiber='B', order_selection=order_selection)

                """ Negative values are just statistical noise around the null flux points,
                removing them would bias the flux level of the sky towards
                higher values

                We proceed in this way:
                -   identify the negative values and make a statistics of their
                    average value
                - if in relevant number, we assume that the their median value
                  corresponds to the noise floor
                - we add the noise floor to the error estimate
                """

                #""" Zero or negative values are identified, flagged and substituted with another value """
                #replacement = 0.01
                observations_B[obs]['null'] = (observations_B[obs]['e2ds'] <= 0.0)
                if (np.sum(observations_B[obs]['null']) > 30):
                    observations_B[obs]['noise_floor'] = np.median(np.abs(
                        observations_B[obs]['e2ds'][observations_B[obs]['null']]))
                else:
                    if observations_B[obs].get('absolute_flux', True):
                        observations_B[obs]['noise_floor'] = 1.0000
                    else:
                        observations_B[obs]['noise_floor'] = 0.00001                #observations_B[obs]['e2ds'][observations_B[obs]['null']] = replacement

                observations_B[obs]['e2ds_err'] = np.sqrt(np.abs(observations_B[obs]['e2ds'])) + observations_B[obs]['noise_floor']

                if 'n_orders' not in observations_B or 'n_pixels' not in observations_B:
                    observations_B['n_orders'] = observations_B[obs]['n_orders']
                    observations_B['n_pixels'] = observations_B[obs]['n_pixels']

                    calib_data_B = get_calib_data(instrument, night_dict[night], archive_dir + night_name, obs,
                                                  fiber='B', order_selection=order_selection)
                has_fiber_B = True
            except:
                pass

        """ Building the base (array of wavelengths) for coadded spectra within the same night """

        observations_A['coadd']['wave'] = np.arange(observations_A['coadd']['wavelength_range'][0],
                                                    observations_A['coadd']['wavelength_range'][1],
                                                    instrument_dict[instrument]['wavelength_step'], dtype=np.double)

        observations_A['coadd']['size'] = np.size(observations_A['coadd']['wave'])
        observations_A['coadd']['step'] = np.ones(observations_A['coadd']['size'], dtype=np.double) * \
                                          instrument_dict[instrument]['wavelength_step']
        print()
        print("  Fixing the observation lists if they are missing")
        if write_transit_list:
            lists_dictionary = _write_transit_list(observations_A,
                                                   lists_dictionary,
                                                   night_dict[night],
                                                   planet_dict)

        print()
        print("  Computing the RV shift outside the transit, and store it to an additional file for quick access ")
        observational_parameters = _get_observational_parameters(observations_A,
                                                                 lists_dictionary,
                                                                 night_dict[night],
                                                                 instrument_dict[instrument],
                                                                 star_dict,
                                                                 planet_dict)
        print()
        for key_name, key_val in observational_parameters['RV_star'].items():
            print(" RV star  {0:s}: {1:f}".format(key_name, key_val))

        print()
        print("  Writing dataset files for night ", night, "  fiber A")
        save_to_cpickle('lists', lists_dictionary, config_in['output'], night)
        save_to_cpickle('input_dataset_fibA', observations_A, config_in['output'], night)
        save_to_cpickle('input_dataset_s1d_fibA', observations_s1d_A, config_in['output'], night)
        save_to_cpickle('calibration_fibA', calib_data_A, config_in['output'], night)
        save_to_cpickle('observational_pams', observational_parameters, config_in['output'], night)

        if has_fiber_B:
            observations_B['coadd'] = {}
            observations_B['coadd']['wave'] = observations_A['coadd']['wave'].copy()
            observations_B['coadd']['size'] = np.size(observations_B['coadd']['wave'])
            observations_B['coadd']['step'] = np.ones(observations_B['coadd']['size'], dtype=np.double) * \
                                              instrument_dict[instrument]['wavelength_step']

            print("  Writing dataset files for night ", night, "  fiber B")
            save_to_cpickle('input_dataset_fibB', observations_B, config_in['output'], night)
            save_to_cpickle('calibration_fibB', calib_data_B, config_in['output'], night)

        """ Running some checks to see if input parameters have been configured properly """
        wavelength_rescaling = instrument_dict[instrument]['wavelength_rescaling']

        pass_wavelength_rescaling = _check_wavelength_rescaling(
            wavelength_rescaling,
            observations_A['coadd']['wavelength_range']
        )
        print()

        if loaded_shared_data: continue

        """ Setting up the base arrays for all the nights"""

        shared_data = _check_coadd_in_shared_data(shared_data, observations_A['coadd']['wavelength_range'])

    if loaded_shared_data: return

    """ Building the base (array of wavelengths) for master-out and coadd spectra
        We do it now to be sure that it will be the same for the whole pipeline
    """
    print("  Creating the shared arrays")
    print()

    shared_data['coadd']['wavelength_range'][0] += 2.0
    shared_data['coadd']['wavelength_range'][1] -= 2.0

    shared_data['coadd']['wave'] = np.arange(shared_data['coadd']['wavelength_range'][0],
                                             shared_data['coadd']['wavelength_range'][1],
                                             config_in['instruments']['shared']['wavelength_step'], dtype=np.double)
    shared_data['coadd']['size'] = np.size(shared_data['coadd']['wave'])
    shared_data['coadd']['step'] = np.ones(shared_data['coadd']['size'], dtype=np.double) * \
                                   config_in['instruments']['shared']['wavelength_step']

    shared_data['binned']['wave'] = np.arange(shared_data['coadd']['wavelength_range'][0],
                                              shared_data['coadd']['wavelength_range'][1],
                                              config_in['instruments']['shared']['wavelength_step']
                                              * config_in['master-out']['binning_factor'], dtype=np.double)
    shared_data['binned']['size'] = np.size(shared_data['binned']['wave'])
    shared_data['binned']['step'] = np.ones(shared_data['binned']['size'], dtype=np.double) * \
                                    config_in['instruments']['shared']['wavelength_step'] * \
                                    config_in['master-out']['binning_factor']

    if config_in['master-out']['wavelength_range'][0] < shared_data['coadd']['wavelength_range'][0] or \
            config_in['master-out']['wavelength_range'][1] > shared_data['coadd']['wavelength_range'][1]:
        warnings.warn("ERROR: Valid master_out wavelength window must be between {0:8.2f} and {1:8.2f}".format(
            shared_data['coadd']['wavelength_range'][0],
            shared_data['coadd']['wavelength_range'][1]))

        pass_wavelength_master_out = False

    shared_data['master-out'] = {}
    shared_data['master-out']['wave'] = np.arange(config_in['master-out']['wavelength_range'][0],
                                                  config_in['master-out']['wavelength_range'][1],
                                                  config_in['master-out']['wavelength_step'], dtype=np.double)

    shared_data['master-out']['size'] = np.size(shared_data['master-out']['wave'])
    shared_data['master-out']['step'] = np.ones(shared_data['master-out']['size'], dtype=np.double) * \
                                        config_in['master-out']['wavelength_step']

    if not (pass_wavelength_master_out and pass_wavelength_rescaling):
        raise ValueError("ERROR: check the previous warnings to see where you are doing it worng")

    print("  COADD wavelength range between {0:8.2f} and {1:8.2f}".format(
        shared_data['coadd']['wavelength_range'][0], shared_data['coadd']['wavelength_range'][1]))
    print("  COADD wavelength step: {0:5.3f}".format(config_in['instruments']['shared']['wavelength_step']))
    print()
    print("Saving shared data")
    save_to_cpickle('shared', shared_data, config_in['output'])
    print()


def _write_transit_list(observations_A, lists_dict, night_dict_key, planet_dict):
    fileout_transit_in_list = open(night_dict_key['in_transit'], 'w')
    fileout_transit_out_list = open(night_dict_key['out_transit'], 'w')
    fileout_transit_full_list = open(night_dict_key['full_transit'], 'w')

    try:
        total_transit_start = np.atleast_1d(night_dict_key['time_of_transit'])[0] - np.atleast_1d(planet_dict['total_transit_duration'])[0] / 2.
        total_transit_end = np.atleast_1d(night_dict_key['time_of_transit'])[0] + np.atleast_1d(planet_dict['total_transit_duration'])[0] / 2.
        full_transit_start = np.atleast_1d(night_dict_key['time_of_transit'])[0] - np.atleast_1d(planet_dict['full_transit_duration'])[0] / 2.
        full_transit_end = np.atleast_1d(night_dict_key['time_of_transit'])[0] + np.atleast_1d(planet_dict['full_transit_duration'])[0] / 2.
    except KeyError:
        total_transit_start = np.atleast_1d(night_dict_key['time_of_transit'])[0] - np.atleast_1d(planet_dict['transit_duration'])[0] / 2.
        total_transit_end = np.atleast_1d(night_dict_key['time_of_transit'])[0] + np.atleast_1d(planet_dict['transit_duration'])[0] / 2.
        full_transit_start = np.atleast_1d(night_dict_key['time_of_transit'])[0] - np.atleast_1d(planet_dict['transit_duration'])[0] / 2.
        full_transit_end = np.atleast_1d(night_dict_key['time_of_transit'])[0] + np.atleast_1d(planet_dict['transit_duration'])[0] / 2.
        print('*** unclear transit duration, ingress/egress observations will be considered full-transit')

    for obs in lists_dict['observations']:

        """ Check if the file should be in transit_in or transit_out list, in case they are not present"""
        #phase_internal = (observations_A[obs]['BJD'] - np.atleast_1d(night_dict_key['time_of_transit'])[0]) / \
        #                 np.atleast_1d(planet_dict['period'])[0]
        #if np.abs(phase_internal) <= planet_dict['transit_duration'][0] / 2. / planet_dict['period'][0]:
        #    fileout_transit_in_list.write('{0:s}\n'.format(obs))
        #else:
        #    fileout_transit_out_list.write('{0:s}\n'.format(obs))

        exptime_seconds = observations_A[obs]['EXPTIME'] / 86400.
        """BJD times have been already corrected to match mid-exposure epochs  """

        if observations_A[obs]['BJD'] + exptime_seconds/2. < total_transit_start \
            or observations_A[obs]['BJD'] - exptime_seconds/2. > total_transit_end:
            fileout_transit_out_list.write('{0:s}\n'.format(obs))
        else:
            fileout_transit_in_list.write('{0:s}\n'.format(obs))
        if observations_A[obs]['BJD'] - exptime_seconds/2. > full_transit_start \
            and observations_A[obs]['BJD'] + exptime_seconds/2. < full_transit_end:
            fileout_transit_full_list.write('{0:s}\n'.format(obs))

    fileout_transit_in_list.close()
    fileout_transit_out_list.close()
    fileout_transit_full_list.close()
    files_list, files_transit_out, files_transit_in, files_transit_full, files_telluric, files_star_telluric = get_filelists(night_dict_key)

    lists_dict['transit_out'] = files_transit_out
    lists_dict['transit_in'] = files_transit_in
    lists_dict['transit_full'] = files_transit_full
    lists_dict['n_transit_out'] = np.size(files_transit_out)
    lists_dict['n_transit_in'] = np.size(files_transit_in)
    lists_dict['n_transit_full'] = np.size(files_transit_full)


    try:
        lists_dict['telluric'] = files_telluric
        lists_dict['n_tellurics'] = np.size(files_telluric)
    except:
        lists_dict['telluric'] = files_transit_out.copy()
        lists_dict['n_tellurics'] = np.size(files_transit_out)

    print(" # observations: ", lists_dict['n_observations'])
    print(" # out-transit obs: ", lists_dict['n_transit_out'])
    print(" # in-transit obs: ", lists_dict['n_transit_in'])
    print(" # full-transit obs: ", lists_dict['n_transit_full'])

    return lists_dict


def _get_observational_parameters(observations_A, lists_dict, night_dict_key, instrument_dict_key, star_dict, planet_dict):
    """
    Compute observational parameters including radial velocity shifts for each observation.
    
    This function calculates various RV shifts needed to transform spectra between
    different reference frames (Observer RF, Solar Barycentric RF, Stellar RF, Planet RF).
    It supports both circular and eccentric planetary orbits.
    
    Parameters
    ----------
    observations_A : dict
        Dictionary containing observation data (BJD, RVC, BERV, AIRMASS, etc.)
        for each observation in fiber A.
    lists_dict : dict
        Dictionary containing lists of observations categorized by transit phase
        (transit_out, transit_in, transit_full, observations, etc.).
    night_dict_key : dict
        Configuration dictionary for the specific night, containing instrument,
        mask, time_of_transit, and other night-specific parameters.
    instrument_dict_key : dict
        Configuration dictionary for the instrument, containing data_archive,
        wavelength_rescaling, refraction parameters, etc.
    star_dict : dict
        Dictionary containing stellar parameters (mass, radius, vsini, etc.)
        and optionally RV_gamma and RV_semiamplitude.
    planet_dict : dict
        Dictionary containing planetary parameters. Required keys:
        - period : [value, error] - Orbital period [days]
        - RV_semiamplitude : [value, error] - Planet's RV semi-amplitude [km/s]
        - reference_time_of_transit : [value, error] - Transit center time [BJD]
        
        For eccentric orbits (orbit != 'circular'), also required:
        - orbit : str - 'circular' or 'eccentric'
        - eccentricity : [value, error] - Orbital eccentricity
        - omega : [value, error] - Argument of periastron [degrees]
        - reference_time : float (optional) - Reference epoch for ephemeris [BJD]
    
    Returns
    -------
    observational_parameters : dict
        Dictionary containing computed parameters for each observation:
        - BJD, mBJD, RVC, AIRMASS, BERV, EXPTIME: raw observation data
        - RV_bjdshift: stellar RV shift due to orbital motion
        - rv_shift_ORF2SRF: RV shift from Observer RF to Stellar RF
        - rv_shift_ORF2BRF: RV shift from Observer RF to Barycentric RF
        - rv_shift_ORF2PRF: RV shift from Observer RF to Planet RF
        - rv_shift_SRF2PRF: RV shift from Stellar RF to Planet RF
        - RV_planet: computed planet radial velocity at observation time
        
        Also includes global parameters:
        - RV_star: dict with stellar RV solution (slope, intercept, systemic)
        - BERV_avg: average barycentric Earth RV
        - orbital_parameters: dict with eccentricity, omega_rad, reference_time
    
    Notes
    -----
    The RV shifts are computed as follows:
    
    1. rv_shift_ORF2SRF: Observer RF -> Stellar RF
       = BERV - (RV_systemic + RV_bjdshift)
       
    2. rv_shift_ORF2PRF: Observer RF -> Planet RF  
       = BERV - RV_systemic - RV_planet
       
    3. rv_shift_SRF2PRF: Stellar RF -> Planet RF
       = RV_bjdshift - RV_planet
    
    For circular orbits: RV_planet = K_planet * sin(2*pi*(t-Tc)/P)
    For eccentric orbits: Full Keplerian solution using compute_planet_RV()
    
    See Also
    --------
    compute_planet_RV : Compute planet RV for circular and eccentric orbits
    """

    observational_parameters = {
        'instrument': night_dict_key['instrument'],
        'mask': night_dict_key['mask'],
        'archive_dir': instrument_dict_key['data_archive'],
        'wavelength_rescaling': instrument_dict_key['wavelength_rescaling'],
        'time_of_transit': np.atleast_1d(night_dict_key['time_of_transit'])[0],
        'refraction_method': instrument_dict_key['refraction']['method'],
        'refraction_fit_order': instrument_dict_key['refraction']['fit_order'],
        'refraction_fit_iters': instrument_dict_key['refraction']['fit_iters'],
        'refraction_fit_sigma': instrument_dict_key['refraction']['fit_sigma'],
        'refraction_knots_spacing': instrument_dict_key['refraction']['knots_spacing'],
        'linear_fit_method': instrument_dict_key['linear_fit_method'],
        'n_orders': observations_A['n_orders'],
        'n_pixels': observations_A['n_pixels'],
        'RV_star': {}
    }

    rv_out = []
    bjd0_out = []

    for obs in lists_dict['transit_out']:
        bjd0_out.extend([observations_A[obs]['BJD'] - observational_parameters['time_of_transit']])
        rv_out.extend([observations_A[obs]['RVC']])

    observational_parameters['RV_star']['slope'], \
    observational_parameters['RV_star']['intercept'], \
    observational_parameters['RV_star']['r_value'], \
    observational_parameters['RV_star']['p_value'], \
    observational_parameters['RV_star']['std_err'] = sci_stats.linregress(bjd0_out, rv_out)

    berv_list = []
    rvc_stack = np.zeros(len(lists_dict['observations']))
    for i, obs in enumerate(lists_dict['observations']):
        berv_list.extend([observations_A[obs]['BERV']])
        observational_parameters[obs] = {
            'BJD': observations_A[obs]['BJD'],
            'mBJD': observations_A[obs]['BJD']-2450000.0000,
            'RVC': observations_A[obs]['RVC'],
            'AIRMASS': observations_A[obs]['AIRMASS'],
            'BERV': observations_A[obs]['BERV'],
            'EXPTIME': observations_A[obs]['EXPTIME']
        }
        rvc_stack[i] = observations_A[obs]['RVC']

    observational_parameters['BERV_avg'] = np.average(berv_list)

    observational_parameters['RV_star']['RV_from_CCF'] = False
    observational_parameters['RV_star']['RV_from_analytical_solution'] = False


    if night_dict_key['use_rv_from_ccf']:
        observational_parameters['RV_star']['RV_from_CCF'] = True
        rvc_systemic = np.average(rvc_stack)
    elif night_dict_key['use_analytical_rvs']:
        observational_parameters['RV_star']['RV_from_analytical_solution'] = True
        rvc_systemic = star_dict['RV_gamma'][0]
        observational_parameters['RV_star']['RV_semiamplitude'] = star_dict['RV_semiamplitude'][0]
    else:
        rvc_systemic = observational_parameters['RV_star']['intercept']

    observational_parameters['RV_star']['RV_systemic'] = rvc_systemic

    # ==========================================================================
    # Extract orbital parameters for planetary RV computation
    # Supports both circular (e=0) and eccentric (e>0) orbits
    # ==========================================================================
    if planet_dict.get('orbit', 'circular') == 'circular':
        # Circular orbit: eccentricity is zero, omega is pi/2 by convention
        eccentricity = 0.0
        omega_rad = np.pi / 2.0
        reference_time = planet_dict['reference_time_of_transit'][0]
    else:
        # Eccentric orbit: read eccentricity and omega from config
        eccentricity = planet_dict['eccentricity'][0]
        omega_rad = planet_dict['omega'][0] * np.pi / 180.0  # Convert degrees to radians
        # Use reference_time if provided, otherwise fall back to transit time
        reference_time = planet_dict.get('reference_time', 
                                         planet_dict['reference_time_of_transit'][0])
    
    # Store orbital parameters in observational_parameters for later use
    observational_parameters['orbital_parameters'] = {
        'eccentricity': eccentricity,
        'omega_rad': omega_rad,
        'reference_time': reference_time,
        'orbit_type': planet_dict.get('orbit', 'circular')
    }

    for obs in lists_dict['observations']:

        if night_dict_key['use_rv_from_ccf']:
            rvc_bjdshift = observations_A[obs]['RVC']
        elif night_dict_key['use_analytical_rvs']:
            rvc_bjdshift = - observational_parameters['RV_star']['RV_semiamplitude'] * np.sin(2 * np.pi * \
                           (observational_parameters[obs]['BJD'] - observational_parameters['time_of_transit'])/planet_dict['period'][0])
        else:
            rvc_bjdshift = observational_parameters['RV_star']['slope'] * \
                           (observational_parameters[obs]['BJD'] - observational_parameters['time_of_transit'])

        observational_parameters[obs]['RV_bjdshift'] = rvc_bjdshift
        observational_parameters[obs]['rv_shift_ORF2SRF'] = observational_parameters[obs]['BERV'] - \
                                                            (rvc_systemic + rvc_bjdshift)

        """ Old definition
        observational_parameters[obs]['rv_shift_ORF2SRF'] = observational_parameters[obs]['BERV'] - \
                                                            (observational_parameters['RV_star']['intercept'] +
                                                             observational_parameters['RV_star']['slope'] *
                                                             (observational_parameters[obs][
                                                                  'BJD'] - time_of_transit))
        """

        """ Slight modification of the RV shift to minimize the rebinning error at the wings of the spectra
        BRF = Solar System Barycentric Reference frame
            rv_shift_ORF2BRF = rv_shift_ORF2SRF_mod + rv_shift_ORF2SRF_res
        """
        observational_parameters[obs]['rv_shift_ORF2BRF'] = \
            observational_parameters[obs]['BERV']

        observational_parameters[obs]['rv_shift_ORF2BRF_mod'] = \
            observational_parameters[obs]['BERV'] - observational_parameters['BERV_avg']

        """ Slight modification of the RV shift to minimize the rebinning error at the wings of the spectra
            rv_shift_ORF2SRF = rv_shift_ORF2SRF_mod + rv_shift_ORF2SRF_res
        """
        observational_parameters[obs]['rv_shift_ORF2SRF_mod'] = \
            observational_parameters[obs]['BERV'] - observational_parameters['BERV_avg'] - rvc_bjdshift

        observational_parameters[obs]['rv_shift_ORF2SRF_res'] = \
            observational_parameters['BERV_avg'] - rvc_systemic

        # ======================================================================
        # RV shift from the Observer RF to the Planet RF
        # ======================================================================
        #STRONG ASSUMPTIONS:
        # - there is only the transiting planet in the system
        # - linear approximation or the orbit near the transit event 
        
        # Computation is performed by:
        # 1. Moving to the Solar System Barycenter (BERV correction)
        # 2. Moving to the Stellar System Barycenter (subtracting systemic RV)
        # 3. Moving onto the planet (subtracting planet's RV)
        #
        # For circular orbits: RV_planet = K_planet * sin(2*pi*t/P)
        # For eccentric orbits: full Keplerian solution is used
        # ======================================================================
        
        # Compute planet's radial velocity at the observation time
        # Using the generalized Keplerian function that handles both circular
        # and eccentric orbits
        rv_planet = compute_planet_RV(
            BJD=observational_parameters[obs]['BJD'],
            time_of_transit=observational_parameters['time_of_transit'],
            period=planet_dict['period'][0],
            K_planet=planet_dict['RV_semiamplitude'][0],
            e0=eccentricity,
            omega0=omega_rad,
            reference_time=reference_time,
            output_unit='km/s'
        )
        
        # Store the planet RV for reference
        observational_parameters[obs]['RV_planet'] = rv_planet
        
        # RV shift from Observer RF to Planet RF
        # BERV: Earth's barycentric velocity (Observer -> Solar Barycenter)
        # rvc_systemic: Star's systemic velocity (Solar Barycenter -> Stellar Barycenter)
        # rv_planet: Planet's velocity (Stellar Barycenter -> Planet RF)
        observational_parameters[obs]['rv_shift_ORF2PRF'] = \
            observational_parameters[obs]['BERV'] \
            - rvc_systemic \
            - rv_planet

        # ======================================================================
        # RV shift from Stellar Rest Frame to Planetary Rest Frame
        # ======================================================================
        # Takes into account the RV of star relatively to the Barycenter
        # rvc_bjdshift: star's orbital motion around system barycenter
        # rv_planet: planet's orbital motion
        observational_parameters[obs]['rv_shift_SRF2PRF'] = \
            + rvc_bjdshift \
            - rv_planet

    observational_parameters['rv_shift_ORF2SRF_res'] = \
        observational_parameters['BERV_avg'] - rvc_systemic

    return observational_parameters


def plot_dataset(config_in, night_input=''):
    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_dict:

        print()
        print("Plotting dataset                           Night: ", night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        instrument = night_dict[night]['instrument']
        wavelength_rescaling = instrument_dict[instrument]['wavelength_rescaling']

        """ Retrieving the observations"""
        input_data = load_from_cpickle('input_dataset_fibA', config_in['output'], night)

        """ Retrieving the calibration data """
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)

        colors, cmap, line_colors = make_color_array(lists, input_data)
        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)

        cbax1 = plt.subplot(gs[:, 1])


        for i, obs in enumerate(lists['observations']):

            rescaling = compute_rescaling(input_data[obs]['wave'], input_data[obs]['e2ds'], wavelength_rescaling)
            for order in range(0,input_data[obs]['n_orders']):


                ax1.plot(input_data[obs]['wave'][order,:], input_data[obs]['e2ds'][order,:]/ rescaling, zorder=i, lw=1,
                        c=line_colors[i], alpha=0.5)
                ax2.plot(input_data[obs]['wave'][order,:], input_data[obs]['e2ds'][order,:] / rescaling, zorder=-i, lw=1,
                        c=line_colors[i], alpha=0.5)

        ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
