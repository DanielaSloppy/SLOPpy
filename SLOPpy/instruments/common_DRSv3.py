import numpy as np
from astropy.io import fits
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.constants import *
from SLOPpy.subroutines.common import *

def DRSv3_map_orders_AB(properties, order_selection):

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


def DRSv3_give_back_selected_orders(properties, fiber, order_selection):

    map_orders_A, map_orders_B = DRSv3_map_orders_AB(properties, order_selection)
    if fiber != 'A':
        return map_orders_B
    else:
        return map_orders_A

def DRSv3_get_calib_data(archive, file_rad, keywords, properties, fiber='A', order_selection=None):

    calib_dict = {}

    selected_orders = DRSv3_give_back_selected_orders(properties, fiber, order_selection)

    map_orders_A, map_orders_B = DRSv3_map_orders_AB(properties, order_selection)

    if fiber=='A':
        calib_dict['fibAB_orders_match'] = map_orders_A - np.min(selected_orders)
        """ The first selected order could not have a match in fiber B, so we need to renumber from the first order of
            the input selection, not from the first order that had a match """
    else:
        calib_dict['fibAB_orders_match'] = map_orders_B - np.min(map_orders_B)
        """ Only the orders with a match in fiber A are read in the first place, so we can safely rescale with respect
            to the number of the first order in the matched list """

    e2ds_fits = fits.open(archive+'/'+file_rad+'_e2ds_'+fiber+'.fits')

    if e2ds_fits[0].header[keywords['header_dpr_catg']] != 'SCIENCE':
        return

    try:
        blaze_file = e2ds_fits[0].header[keywords['header_blaze']]
        blaze_fits = fits.open(archive + '/' + blaze_file)
    except:
        blaze_file = e2ds_fits[0].header[keywords['header_blaze']].replace(':', '-')
        blaze_fits = fits.open(archive + '/' + blaze_file)

    # getting blaze file
    calib_dict['blaze'] = blaze_fits[0].data[selected_orders, :]

    # getting lamp file
    try:
        lamp_fits = fits.open(archive + '/' + blaze_file[:29] + '_lamp_' + fiber + '.fits')
        calib_dict['lamp'] = lamp_fits[0].data[selected_orders, :]
        lamp_fits.close()
    except:
        print("lamp files not available, sky correction will not be performed")

    calib_dict['n_pixels'] = blaze_fits[0].header['NAXIS1']
    calib_dict['n_orders'] = len(selected_orders)

    blaze_fits.close()

    return calib_dict


def DRSv3_get_input_data(archive, file_rad, keywords, properties, mask, fiber='A', skip_ccf=None, skip_s1d=True, order_selection=None):

    input_dict = {'mask': mask, 'header':{}}
    input_s1d = {'header':{}}

    selected_orders = DRSv3_give_back_selected_orders(properties, fiber, order_selection)

    if mask is None:
        skip_ccf = True

    e2ds_fits = fits.open(archive+'/'+file_rad+'_e2ds_'+fiber+'.fits')
    input_dict['header']['e2ds'] = e2ds_fits[0].header

    input_dict['n_pixels'] = e2ds_fits[0].header['NAXIS1']
    input_dict['n_orders'] = len(selected_orders)

    input_dict['DPR_CATG'] = e2ds_fits[0].header[keywords['header_dpr_catg']]

    if input_dict['DPR_CATG'] != 'SCIENCE':
        return
    input_dict['DPR_TYPE'] = e2ds_fits[0].header[keywords['header_dpr_type']]

    if not skip_s1d:
        s1d_fits = fits.open(archive + '/' + file_rad + '_s1d_'+fiber+'.fits')
        input_dict['header']['s1d'] = s1d_fits[0].header
        input_s1d['header']['s1d'] = s1d_fits[0].header

        temp_wave, temp_step = DRSv3_get_s1d_wave(s1d_fits)
        sel_wave = (temp_wave >= 3879.99990) & (temp_wave <= 6900.0001)

        input_s1d['flux'] = s1d_fits[0].data[sel_wave]
        input_s1d['wave'] = temp_wave[sel_wave]
        input_s1d['step'] = temp_step[sel_wave]
        input_s1d['size'] = np.size(input_s1d['wave'])

        s1d_fits.close()

    if not skip_ccf:
        ccf_fits = fits.open(archive+'/'+file_rad+'_ccf_'+mask+'_'+fiber+'.fits')
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

    if properties['time_standard'] == 'UTC':
        input_dict['BJD']+= difference_utc2tdb(input_dict['MJD']+2400000.5)


    input_dict['LST'] = e2ds_fits[0].header['LST']

    try:
        input_dict['AIRMASS'] = e2ds_fits[0].header['AIRMASS']
    except:
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

        try:
            input_dict['AIRMASS'] = e2ds_fits[0].header['AIRMASS']
        except:
            input_dict['AIRMASS'] = (e2ds_fits[0].header[keywords['airmass_alt_start']]
										+ e2ds_fits[0].header[keywords['airmass_alt_end']])/2.

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

    except:
        print('Keyword error in prepare_dataset - check the FITS header of your files')
        quit()
        pass

    input_dict['BLAZE_file'] = e2ds_fits[0].header[keywords['header_blaze']]
    input_dict['CCD_SIGDET'] = e2ds_fits[0].header[keywords['header_ccd']]
    input_dict['CCD_GAIN'] = e2ds_fits[0].header[keywords['header_conad']]

    # getting data
    input_dict['e2ds'] = e2ds_fits[0].data[selected_orders, :]
    input_dict['e2ds_err'] = np.sqrt(np.abs(input_dict['e2ds']))

    temp_wave, temp_step = DRSv3_get_e2ds_wave(e2ds_fits, keywords['header_deg_ll'], keywords['header_coeff_ll'])

    input_dict['wave'] = temp_wave[selected_orders, :]
    input_dict['step'] = temp_step[selected_orders, :]
    input_dict['orders'] = len(selected_orders)

    input_dict['wave_size'] = e2ds_fits[0].header['NAXIS1']

    e2ds_fits.close()
    if not skip_ccf:
        ccf_fits.close()
        
    input_dict['e2ds']=cleanSpectrumMatrix(input_dict['wave'],input_dict['e2ds'])
    return input_dict,input_s1d


def DRSv3_get_s1d_wave(s1d_fits):
    return np.arange(0, s1d_fits[0].header['NAXIS1'], 1.)*s1d_fits[0].header['CDELT1'] + s1d_fits[0].header['CRVAL1'], \
           np.ones(s1d_fits[0].header['NAXIS1'])*s1d_fits[0].header['CDELT1']


def DRSv3_get_e2ds_wave(e2ds_fits, header_deg_ll, header_coeff_ll, order=None):

    e2ds_o = e2ds_fits[0].header['NAXIS2']
    e2ds_w = e2ds_fits[0].header['NAXIS1']

    e2ds_wave = np.zeros([e2ds_o, e2ds_w], dtype=np.double)
    e2ds_step = np.zeros([e2ds_o, e2ds_w], dtype=np.double)

    d = e2ds_fits[0].header[header_deg_ll]
    x = np.arange(0, e2ds_w, 1.)

    for n in range(0, e2ds_o):
        for i in range(d, -1, -1):
            a_sel = i + n*(1+d)
            a_coeff = e2ds_fits[0].header[header_coeff_ll+repr(a_sel)]
            if i == d:
                y_w = a_coeff
                y_s = i*a_coeff
            else:
                y_w = y_w*x + a_coeff
                if i > 0: y_s = y_s*x + i*a_coeff
        e2ds_wave[n, :] = y_w
        e2ds_step[n, :] = y_s

    if order is None:
        return e2ds_wave, e2ds_step
    else:
        return e2ds_wave[order, :], e2ds_step[order, :]
