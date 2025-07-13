import numpy as np
from astropy.io import fits
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.constants import *
from SLOPpy.subroutines.common import *

def ESPRESSO_map_orders_AB(properties, order_selection):

    n_orders_A = properties['n_orders_A']
    map_orders_A_full = np.arange(0, n_orders_A, dtype=np.int16)
    map_orders_BtoA = np.asarray(properties['orders_BtoA'])

    map_orders_A = []
    map_orders_B = []

    #if order_selection == 'red_ccd':
    #    working_A_full = map_orders_A_full[properties['red_ccd']]
    #    working_BtoA = map_orders_BtoA[properties['red_ccd']]
    if order_selection in ['red_ccd', 'blue_ccd', 'full_ccd']:
        working_A_full = map_orders_A_full[properties[order_selection]]
        working_BtoA = map_orders_BtoA[properties[order_selection]]

    elif order_selection:
        working_A_full = map_orders_A_full[order_selection]
        working_BtoA = map_orders_BtoA[order_selection]
    else:
        working_A_full = map_orders_A_full
        working_BtoA = map_orders_BtoA

    for order_A, order_B in zip(working_A_full, working_BtoA):
        if order_B < 0: continue
        map_orders_A.extend([order_A])
        map_orders_B.extend([order_B])

    return map_orders_A, map_orders_B


def ESPRESSO_give_back_selected_orders(properties, fiber, order_selection):

    map_orders_A, map_orders_B = ESPRESSO_map_orders_AB(properties, order_selection)
    if fiber != 'A':
        return map_orders_B
    else:
        return map_orders_A

def ESPRESSO_get_calib_data(archive, file_rad, night_dict, fiber='A', order_selection=None):

    calib_dict = {}

    keywords, properties = ESPRESSO_get_instrument_keywords(night_dict)

    selected_orders = ESPRESSO_give_back_selected_orders(properties, fiber, order_selection)

    map_orders_A, map_orders_B = ESPRESSO_map_orders_AB(properties, order_selection)

    if fiber=='A':
        calib_dict['fibAB_orders_match'] = map_orders_A - np.min(selected_orders)
        """ The first selected order could not have a match in fiber B, so we need to renumber from the first order of
            the input selection, not from the first order that had a match """
    else:
        calib_dict['fibAB_orders_match'] = map_orders_B - np.min(map_orders_B)
        """ Only the orders with a match in fiber A are read in the first place, so we can safely rescale with respect
            to the number of the first order in the matched list """

    if properties['use_ESO_deblazed']:
        if properties['use_ESO_telluric_correction']:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_TELL_CORR_'+fiber+'.fits')
        elif properties['use_ESO_sky_correction']:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_SKYSUB_'+fiber+'.fits')
        else:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_'+fiber+'.fits')

        blaze = np.ones_like(e2ds_fits[1].data[properties['orders_regroup'], :])
        e2ds_fits.close()
    else:
        e2ds_fits_deblazed = fits.open(archive+'/'+file_rad+'_S2D_'+fiber+'.fits')
        e2ds_fits= fits.open(archive+'/'+file_rad+'_S2D_BLAZE_'+fiber+'.fits')

        blaze_full = e2ds_fits[1].data / e2ds_fits_deblazed[1].data
        blaze = blaze_full[properties['orders_regroup'], :]

        e2ds_fits_deblazed.close()
        e2ds_fits.close()

    calib_dict['blaze'] = blaze[selected_orders, :]

    # getting lamp file
    #try:
    #    lamp_fits = fits.open(archive + '/' + blaze_file[:29] + '_lamp_' + fiber + '.fits')
    #    calib_dict['lamp'] = lamp_fits[0].data[selected_orders, :]
    #    lamp_fits.close()
    #except:
    #    print("lamp files not available, sky correction cannot not be performed")

    calib_dict['n_pixels'] = np.shape(blaze)[1]
    calib_dict['n_orders'] = len(selected_orders)

    return calib_dict


def ESPRESSO_get_input_data(archive, file_rad, night_dict, fiber='A', skip_ccf=None, skip_s1d=True, order_selection=None):

    input_dict = {'header':{}}
    input_s1d = {'header':{}}

    keywords, properties = ESPRESSO_get_instrument_keywords(night_dict)

    selected_orders = ESPRESSO_give_back_selected_orders(properties, fiber, order_selection)

    if properties['use_ESO_deblazed']:
        if properties['use_ESO_telluric_correction']:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_TELL_CORR_'+fiber+'.fits')
        elif properties['use_ESO_sky_correction']:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_SKYSUB_'+fiber+'.fits')
        else:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_'+fiber+'.fits')
    else:
        if properties['use_ESO_telluric_correction']:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_BLAZE_TELL_CORR_'+fiber+'.fits')
        else:
            e2ds_fits = fits.open(archive+'/'+file_rad+'_S2D_BLAZE_'+fiber+'.fits')

    input_dict['header']['e2ds'] = e2ds_fits[0].header

    input_dict['n_pixels'] = e2ds_fits[1].header['NAXIS1']
    input_dict['n_orders'] = len(selected_orders)

    input_dict['DPR_CATG'] = e2ds_fits[0].header[keywords['header_dpr_catg']]
    if e2ds_fits[0].header[keywords['header_dpr_catg']] != True:
        return

    if not skip_s1d:
        if properties['use_ESO_telluric_correction']:
            s1d_fits = fits.open(archive + '/' + file_rad + '_S1D_TELL_CORR_'+fiber+'.fits')
        elif properties['use_ESO_sky_correction']:
            s1d_fits = fits.open(archive + '/' + file_rad + '_S1D_SKYSUB_'+fiber+'.fits')
        else:
            s1d_fits = fits.open(archive + '/' + file_rad + '_S1D_'+fiber+'.fits')

        input_dict['header']['s1d'] = s1d_fits[0].header
        input_s1d['header']['s1d'] = s1d_fits[0].header

        temp_wave = s1d_fits[1].data['wavelength_air']
        temp_step = np.zeros_like(temp_wave)
        temp_step[1:] = np.diff(temp_wave)
        temp_step[0] = temp_step[1]
        sel_wave = (temp_wave >= 3799.9999) & (temp_wave <= 7850.0001)

        input_s1d['flux'] = s1d_fits[1].data['flux'][sel_wave]
        input_s1d['flux_err'] = s1d_fits[1].data['error'][sel_wave]
        input_s1d['wave'] = temp_wave[sel_wave]
        input_s1d['step'] = temp_step[sel_wave]
        input_s1d['size'] = np.size(input_s1d['wave'])

        s1d_fits.close()

    if not skip_ccf:
        if properties['use_ESO_telluric_correction']:
            ccf_fits = fits.open(archive+'/'+file_rad+'_CCF_TELL_CORR_'+fiber+'.fits')
        elif properties['use_ESO_sky_correction']:
            ccf_fits = fits.open(archive+'/'+file_rad+'_CCF_SKYSUB_'+fiber+'.fits')
        else:
            ccf_fits = fits.open(archive+'/'+file_rad+'_CCF_'+fiber+'.fits')
        input_dict['RVC'] = ccf_fits[0].header[keywords['header_rvc']]
        input_dict['header']['ccf'] = ccf_fits[0].header



    input_dict['BERV'] = e2ds_fits[0].header[keywords['header_berv']]
    input_dict['EXPTIME'] = e2ds_fits[0].header['EXPTIME']

    if properties['time_stamp'] == 'start_exposure':
        input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']] + input_dict['EXPTIME']/86400.
        input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']] + input_dict['EXPTIME']/86400.
    elif properties['time_stamp'] == 'mid_exposure':
        input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']]
        input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']]
    elif properties['time_stamp'] == 'end_exposure':
        input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']] - input_dict['EXPTIME']/86400.
        input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']] - input_dict['EXPTIME']/86400.
    else:
        print('*** please specify the relationship between epoch and exposure time - assuming mid-exposure epochs')
        input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']]
        input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']]

    input_dict['LST'] = e2ds_fits[0].header['LST']

    input_dict['AIRMASS'] = (e2ds_fits[0].header[keywords['airmass_alt_start']]
                                + e2ds_fits[0].header[keywords['airmass_alt_end']])/2.

    try:
        input_dict['AIRM_START'] = e2ds_fits[0].header[keywords['airmass_alt_start']]
        input_dict['AIRM_END'] = e2ds_fits[0].header[keywords['airmass_alt_end']]
    except:
        input_dict['AIRM_START'] = input_dict['AIRMASS'] -0.05
        input_dict['AIRM_END'] = max(1.00, input_dict['AIRMASS'] +0.05)

    input_dict['UTC'] = (input_dict['MJD'] - int(input_dict['MJD'])) * 86400.
    input_dict['HUMIDITY'] = e2ds_fits[0].header[keywords['humidity']]
    input_dict['PRESSURE'] = e2ds_fits[0].header[keywords['pressure']]
    input_dict['TEMPERATURE_EN'] = e2ds_fits[0].header[keywords['temperature_env']]
    input_dict['TEMPERATURE_M1'] = e2ds_fits[0].header[keywords['temperature_m1']]
    input_dict['ELEVATION'] = np.arcsin(1./input_dict['AIRMASS']) * (180./np.pi)

    input_dict['GEOELEV'] = properties['geoelev']
    input_dict['GEOLONG'] = properties['longitude']
    input_dict['GEOLAT'] = properties['latitude']

    input_dict['molecfit'] = properties['molecfit']

    try:
        try:
            input_dict['RA'] = e2ds_fits[0].header['RA-DEG']
            input_dict['DEC'] = e2ds_fits[0].header['DEC-DEG']
        except:
            input_dict['RA'] = e2ds_fits[0].header['RA-RAD'] * 180.00 / np.pi
            input_dict['DEC'] = e2ds_fits[0].header['DEC-RAD'] * 180.00 / np.pi  # weird choice of using DEC in hours
    except:
            input_dict['RA'] = e2ds_fits[0].header['RA']
            input_dict['DEC'] = e2ds_fits[0].header['DEC']

    try:
        input_dict['BERV'] = e2ds_fits[0].header[keywords['header_berv']]
        input_dict['EXPTIME'] = e2ds_fits[0].header['EXPTIME']

        if properties['time_stamp'] == 'start_exposure':
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']] + input_dict['EXPTIME']/86400.
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']] + input_dict['EXPTIME']/86400.
        elif properties['time_stamp'] == 'mid_exposure':
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']]
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']]
        elif properties['time_stamp'] == 'end_exposure':
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']] - input_dict['EXPTIME']/86400.
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']] - input_dict['EXPTIME']/86400.
        else:
            print('*** please specify the relationship between epoch and exposure time - assuming mid-exposure epochs')
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']]
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']]

        if properties['time_standard'] == 'UTC':
            input_dict['BJD']+= difference_utc2tdb(input_dict['MJD']+2400000.5)

        input_dict['LST'] = e2ds_fits[0].header['LST']


        input_dict['AIRMASS'] = (e2ds_fits[0].header[keywords['airmass_alt_start']]
                                    + e2ds_fits[0].header[keywords['airmass_alt_end']])/2.

        input_dict['UTC'] = e2ds_fits[0].header[keywords['header_utc']]
        #input_dict['UTC'] = (input_dict['MJD'] - int(input_dict['MJD'])) * 86400.
        input_dict['HUMIDITY'] = e2ds_fits[0].header[keywords['humidity']]
        input_dict['PRESSURE'] = e2ds_fits[0].header[keywords['pressure']]
        input_dict['TEMPERATURE_EN'] = e2ds_fits[0].header[keywords['temperature_env']]
        input_dict['TEMPERATURE_M1'] = e2ds_fits[0].header[keywords['temperature_m1']]
        input_dict['ELEVATION'] = np.arcsin(1./input_dict['AIRMASS']) * (180./np.pi)

        input_dict['GEOELEV'] = properties['geoelev']
        input_dict['GEOLONG'] = properties['longitude']
        input_dict['GEOLAT'] = properties['latitude']

        input_dict['molecfit'] = properties['molecfit']

        input_dict['RA'] = e2ds_fits[0].header['RA']
        input_dict['DEC'] = e2ds_fits[0].header['DEC']

        """""
        try:
            try:
                input_dict['RA'] = e2ds_fits[0].header['RA-DEG']
                input_dict['DEC'] = e2ds_fits[0].header['DEC-DEG']
            except:
                input_dict['RA'] = e2ds_fits[0].header['RA-RAD'] * 180.00 / np.pi
                input_dict['DEC'] = e2ds_fits[0].header['DEC-RAD'] * 180.00 / np.pi  # weird choice of using DEC in hours
        except:
                input_dict['RA'] = e2ds_fits[0].header['RA']
                input_dict['DEC'] = e2ds_fits[0].header['DEC']
        """""
    except:
        print('Keyword error in prepare_dataset - check the FITS header of your files')
        quit()
        pass

    # getting data

    input_dict['e2ds'] = e2ds_fits[1].data[properties['orders_regroup'], :][selected_orders, :]
    input_dict['e2ds_err'] = e2ds_fits[2].data[properties['orders_regroup'], :][selected_orders, :]

    if properties['apply_ESO_telluric_correction']:
        ssub_fits = fits.open(archive+'/'+file_rad+'_S2D_SKYSUB_'+fiber+'.fits')
        tell_fits = fits.open(archive+'/'+file_rad+'_S2D_TELL_CORR_'+fiber+'.fits')
        ssub_spec = ssub_fits[1].data[properties['orders_regroup'], :][selected_orders, :]
        tell_spec = tell_fits[1].data[properties['orders_regroup'], :][selected_orders, :]

        ssub_spec[ssub_spec == 0] = 1.0
        tell_spec[tell_spec == 0] = 1.0

        tell_spectrum_selected = tell_spec / ssub_spec
        ssub_fits.close()
        tell_fits.close()

        tell_spectrum_selected[tell_spectrum_selected == 0] = 1.0
        input_dict['e2ds'] /= tell_spectrum_selected
        input_dict['e2ds_err'] /= tell_spectrum_selected

    wave_with_berv =  e2ds_fits[5].data[properties['orders_regroup'], :][selected_orders, :]
    step_with_berv =  e2ds_fits[7].data[properties['orders_regroup'], :][selected_orders, :]

    #input_dict['wave'], input_dict['step'] = shift_wavelength_to_rest(wave_with_berv, step_with_berv, -input_dict['BERV'] )
    input_dict['wave'], input_dict['step'] = shift_wavelength(wave_with_berv, step_with_berv, -input_dict['BERV'] )
    input_dict['wave_size'] = e2ds_fits[1].header['NAXIS1']

    e2ds_fits.close()


    if not skip_ccf:
        ccf_fits.close()

    return input_dict,input_s1d



def ESPRESSO_get_instrument_keywords(night_dict):

    if night_dict['telescope'] == 'UT1':
        telescope = '1'
    elif night_dict['telescope'] == 'UT2':
        telescope = '2'
    elif night_dict['telescope'] == 'UT3':
        telescope = '3'
    elif night_dict['telescope'] == 'UT4':
        telescope = '4'
    else:
        telescope = ''

    keywords = {
        'header_rvc': 'HIERARCH ESO QC CCF RV',
        'header_berv': 'HIERARCH ESO QC BERV',
        'header_bjd': 'HIERARCH ESO QC BJD',
        'header_mjd': 'MJD-OBS', # MJD in days. This parameter is required for the retrieval of GDAS data
        'header_utc': 'UTC',

        'header_dpr_catg': 'HIERARCH ESO PRO SCIENCE',

        'airmass_alt_start': 'HIERARCH ESO TEL'+telescope+' AIRM START',
        'airmass_alt_end': 'HIERARCH ESO TEL'+telescope+' AIRM END',
        ## Telescope altitude is computed using the middle values obtained from airmass

        'humidity':'HIERARCH ESO TEL'+telescope+' AMBI RHUM', # Relative humidity in % for GEOELEV.
        #'pressure_start' : 'HIERARCH ESO TEL AMBI PRES START',
        #'pressure_end': 'HIERARCH ESO TEL AMBI PRES END',
        'pressure':'HIERARCH ESO TEL'+telescope+' AMBI PRES END',
        'temperature_env': 'HIERARCH ESO TEL'+telescope+' AMBI TEMP', #Ambient temperature in C for GEOELEV
        'temperature_m1': 'HIERARCH ESO TEL'+telescope+' TH M1 TEMP', # Temperature of primary mirror M1 in C (for emission spectra only)

    }

    properties = {
        # DRS-specific keywords
        'time_stamp': 'mid_exposure',
        'time_standard': 'BJD_TDB',

        # Observatory-specific keywords
        'geoelev': 2648.0, # meters
        'longitude' : -70.4051, # Tel geo longitude (+=East) (deg)
        'latitude' : -24.6276,  # Tel geo latitute (+=North) (deg)
        # Instrument-specific keyword
        'n_orders_A': 85,
        'n_orders_B': 85,
        'orders_BtoA':
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                80, 81, 82, 83, 84],
            # after many experiments, I found out the easiest and more robust way to define
            # the order correspondence between fiber A anf B is just to write it down
        'red_ccd':
            [                               47, 48, 49,
                50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                80, 81, 82, 83, 84],
        'blue_ccd':
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
             10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
             30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
             40, 41, 42, 43, 44, 45, 46],
        'full_ccd':
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
             10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
             30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
             40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
             50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
             60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
             70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
             80, 81, 82, 83, 84],
        # The following are the input values used by Molecfit, taken from Allart+2017
        # for convenience, all the default values are listed here instead of being scattered into the code
        'molecfit': {
            'default_wstep': 0.01000, # default wavelength step size for the input stellar spectra
            'molecules': ['H2O', 'O2'],
            'ftol': 1e-9,
            'xtol': 1e-9,

            'cont_const': 1.0, # a0,  This value differs from Allart+2017 since we are using normalized spectra
            'cont_n': 3, # n_cont, Degree of coefficients for continuum fit
            'wlc_n': 2, # n_lambda, Polynomial degree of the refined wavelength solution
            'wlc_const': 0.0, # b0, Initial constant term for wavelength correction (shift relative to half  wavelength range)
            'res_gauss': 4.8, # omega_Gaussian,  Initial value for FWHM of Gaussian in pixels
            'kernfac': 15, #kernel_size, Size of Gaussian/Lorentzian/Voigtian kernel in FWHM

            'slitwidth': 1.00, # in arcseconds
            'pixelscale': 0.16,
        }
    }

    if night_dict['spectral_selection'] == 'odd':
        properties['orders_regroup'] = [
            1, 3, 5, 7, 9, 11, 13, 15, 17, 19,
            21, 23, 25, 27, 29, 31, 33, 35, 37, 39,
            41, 43, 45, 47, 49, 51, 53, 55, 57, 59,
            61, 63, 65, 67, 69, 71, 73, 75, 77, 79,
            81, 83, 85, 87, 89, 91, 93, 95, 97, 99,
            101, 103, 105, 107, 109, 111, 113, 115, 117, 119,
            121, 123, 125, 127, 129, 131, 133, 135, 137, 139,
            141, 143, 145, 147, 149, 151, 153, 155, 157, 159,
            161, 163, 165, 167, 169]

    if night_dict['spectral_selection'] == 'even':
        properties['orders_regroup'] = [
            0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
            20, 22, 24, 26, 28, 30, 32, 34, 36, 38,
            40, 42, 44, 46, 48, 50, 52, 54, 56, 58,
            60, 62, 64, 66, 68, 70, 72, 74, 76, 78,
            80, 82, 84, 86, 88, 90, 92, 94, 96, 98,
            100, 102, 104, 106, 108, 110, 112, 114, 116, 118,
            120, 122, 124, 126, 128, 130, 132, 134, 136, 138,
            140, 142, 144, 146, 148, 150, 152, 154, 156, 158,
            160, 162, 164, 166, 168]

    properties['use_ESO_telluric_correction'] = night_dict.get('use_ESO_telluric_correction', True)
    properties['use_ESO_sky_correction'] = night_dict.get('use_ESO_sky_correction', True)
    properties['use_ESO_deblazed'] = night_dict.get('use_ESO_deblazed', True)
    properties['apply_ESO_telluric_correction'] = night_dict.get('apply_ESO_telluric_correction', False)

    if properties['apply_ESO_telluric_correction'] and properties['use_ESO_sky_correction']:
        print('*** WARNING: DAS telluric correction is performed on spectra with pre-applied sky correction.' \
              '             Please check the DAS documentation.')
        print('    You shoudl use the keyword use_ESO_telluric_correction to use directly ESO telluric correction')
        quit()

    if properties['use_ESO_telluric_correction'] and not properties['use_ESO_sky_correction']:
        print('*** WARNING: pre-applied telluric correction comes with pre-applied sky correction.' \
              '             a DAS file with only telluric correction is not available.' \
              '             Please check the DAS documentation.')
        print('    You can use the keyword apply_ESO_telluric_correction to apply the ESO telluric correction '\
              '    without using the pre-applied sky correction.')
        quit()

    if properties['use_ESO_sky_correction'] and not properties['use_ESO_deblazed']:
        print('*** WARNING: pre-applied sky correction is available only on deblazed files.' \
              '             As the computation of the relative efficiency between the two fibers is not yet implemented,' \
              '             at the moment is not possible to perform sky correction through the usual SLOPpy task.')

    return keywords, properties

