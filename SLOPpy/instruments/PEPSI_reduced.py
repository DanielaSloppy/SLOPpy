import numpy as np
from astropy.io import fits
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines import constants
from SLOPpy.subroutines.common import *

from astropy.coordinates import SkyCoord
from astropy import units as u

GLOBAL_WAVE_REF = None

def PEPSI_get_instrument_keywords():

    properties = {
        # DRS-specific keywords
        'time_stamp': 'mid_exposure',
        'time_standard': 'TDB',

        # Observatory-specific keywords
        'geoelev': 3221.0, # meters
        'longitude' : -7.325938, # Tel geo longitude (+=East) (deg)
        'latitude' : 32.7013083,  # Tel geo latitute (+=North) (deg)

        # Instrument-specific keyword
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
    return properties

def PEPSI_get_calib_data(archive, file_rad, fiber='A', order_selection=None):
    """ There are no calibration files from PEPSI, so this subroutine will 
        provide the dictionary required by SLOPpy to properly work
        "fiber", "order_selection" variables are kept for consistency with the
        main code, but they do not have applicability 
    """
    calib_dict = {}

    # file from which to extract the relevant keywords - it could be a science
    # frame as well
    pepsi_fits = fits.open(archive+'/'+file_rad)
    data_fits = pepsi_fits[1].data['Fun']

    # Blaze file - if only blaze-corrected files are available, set it equal to 1.
    calib_dict['n_pixels'] = len(data_fits)
    calib_dict['n_orders'] = 1

    calib_dict['blaze'] = np.ones((calib_dict['n_orders'], calib_dict['n_pixels']))


    return calib_dict


def PEPSI_get_input_data(archive, file_rad, mask, fiber='A', skip_ccf=None, skip_s1d=True, order_selection=None):
    """_summary_

    Returns:
        _type_: _description_
    """

    """ PEPSI delivers calibrated, rebinned spectra only in a specific format
        so many entry in the dictionary will be empty.
        Given the simplicity of the format, this subroutines can be used as
        template for other instruments
        "fiber", "skip_ccf", "skip_s1d", "order_selection" variables are kept
        for consistency with the main code, but they do not have applicability
    """

    input_dict = {'mask': mask, 'header':{}}
    input_s1d = {'header':{}}

    properties = PEPSI_get_instrument_keywords()


    pepsi_fits = fits.open(archive+'/'+file_rad)
    input_dict['header']['e2ds'] = pepsi_fits[0].header


    #Arg = file_fits[1].data['Arg']
    #Fun = file_fits[1].data['Fun']
    #Var = file_fits[1].data['Var']
    #Mask = file_fits[1].data['Mask']

    input_dict['n_pixels'] = pepsi_fits[1].header['NAXIS2']
    input_dict['n_orders'] = 1

    # Not sure if these keywords are required anywhere
    #input_dict['DPR_CATG'] = 'SCIENCE'
    #input_dict['DPR_TYPE'] = ?? 

    #input_dict['BERV'] = pepsi_fits[0].header['SSBVEL'] / 1000. # in km/s
    #input_dict['RVC'] = pepsi_fits[0].header['RADVEL'] / 1000.
    # RV of the star, it must be provided but it can be bypassed
    input_dict['BERV'] = 0.0
    input_dict['RVC'] = 0.0
    

    input_dict['EXPTIME'] = pepsi_fits[0].header['EXPTIME']

    # BJD provided at midexposure ,no need to check for it
    input_dict['BJD'] = pepsi_fits[0].header['JD-TDB']
    input_dict['MJD'] = pepsi_fits[0].header['JD-OBS'] - constants.MJD

    input_dict['AIRMASS'] = pepsi_fits[0].header['AIRMASS']
    input_dict['AIRM_START'] = pepsi_fits[0].header['AIRMASS']
    input_dict['AIRM_END'] = pepsi_fits[0].header['AIRMASS']

    input_dict['UTC'] = (input_dict['MJD'] - int(input_dict['MJD'])) * 86400.
    input_dict['HUMIDITY'] = pepsi_fits[0].header['LBTH'] # Relative humidity in % for GEOELEV.
    input_dict['PRESSURE'] = pepsi_fits[0].header['LBTP']
    input_dict['TEMPERATURE_EN'] = pepsi_fits[0].header['LBTT']  #Ambient temperature in C for GEOELEV
    input_dict['TEMPERATURE_M1'] = pepsi_fits[0].header['LBTT']  #Temperature of primary mirror M1 in C (for emission spectra only)
    input_dict['ELEVATION'] = np.arcsin(1./input_dict['AIRMASS']) * (180./np.pi)

    input_dict['GEOELEV'] = properties['geoelev']
    input_dict['GEOLONG'] = properties['longitude']
    input_dict['GEOLAT'] = properties['latitude']

    input_dict['molecfit'] = properties['molecfit']

    skycoords = c = SkyCoord(pepsi_fits[0].header['RA'], pepsi_fits[0].header['DEC'], unit=(u.hourangle, u.deg))
    input_dict['RA'] = c.ra.degree
    input_dict['DEC'] = c.dec.degree

    # Not sure if thes evalies are required, it may be possible
    #input_dict['BLAZE_file'] = None
    #input_dict['CCD_SIGDET'] = None
    #input_dict['CCD_GAIN'] = None

    # getting data

    """
    NOTE: PEPSI provides 1D spectra in the stellar reference frame,
    but SLOPpy requires 2D spectra (order by order) in the observer reference
    frame. Thus:
    1) we shift the wavelength from the stellar to the observer reference frame
    2) we transform the array into (1,n_pixels) shaped arrays

    An empty array as required by SLOPpy would have the shape
    np.empty([n_orders, n_pixels])
    """

    global GLOBAL_WAVE_REF
    wave_stellar =  pepsi_fits[1].data['Arg']
    #rvshift = pepsi_fits[0].header['SSTVEL'] / 1000.
    if GLOBAL_WAVE_REF is None:
        GLOBAL_WAVE_REF = wave_stellar.copy()
        
    wave_ref = GLOBAL_WAVE_REF
    flux = pepsi_fits[1].data['Fun']
    err = pepsi_fits[1].data['Var']
    flux_i = np.interp(wave_ref, wave_stellar,flux)
    err_i = np.interp(wave_ref, wave_stellar,err)

    #input_dict['wave_size'] = pepsi_fits[1].header['NAXIS2']
    #input_dict['wave'] = np.reshape(shift_wavelength_array(wave_stellar, rvshift), (1, input_dict['wave_size']))
    #input_dict['wave'] = np.reshape(wave_stellar, (1, input_dict['wave_size']))
    #input_dict['e2ds'] = np.reshape(pepsi_fits[1].data['Fun'], (1, input_dict['wave_size']))
    #input_dict['e2ds_err'] = np.reshape(pepsi_fits[1].data['Var'], (1, input_dict['wave_size']))
    
    input_dict['wave'] = wave_ref.reshape(1,-1)
    input_dict['wave_size'] = len(wave_ref)
    input_dict['e2ds'] = flux_i.reshape(1,-1)
    input_dict['e2ds_err'] = err_i.reshape(1,-1)

    """ PEPSI spectra are normalized to unity, but SLOPpy is expecting to have spectra in absolute counts
        absolute counts mean that a larger step size will have a larger number of counts given the same
        flux density at a specific wavelength.
        PEPSI spectra have been resampled on a non-linear scale and than normalized, so not taking into account
        the bin (or step) size would introduce a deformation during the rebinning phase
        After several tests, I am forced to introduce a flag to force non-preservation of flux at every
        rebinning step across the code

    """
    input_dict['absolute_flux'] = False

    input_dict['step'] = np.zeros_like(input_dict['wave'])
    input_dict['step'][0,1:-1] = (input_dict['wave'][0,2:] - input_dict['wave'][0,:-2])/2.
    input_dict['step'][0,0] = input_dict['step'][0,1]
    input_dict['step'][0,-1] = input_dict['step'][0,-2]

    # order selection is always equal the the first - and unique - order
    input_dict['orders'] = [0]

    input_dict['absolute_flux'] = False


    pepsi_fits.close()

    return input_dict,input_s1d
