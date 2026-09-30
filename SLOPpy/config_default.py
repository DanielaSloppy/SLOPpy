config_default = {
    'molecfit': {
        'installation_path': '/usr/local/eso/bin/',
        'include_telluric': './additional_files/include_regions_ORF.dat',
        'include_stellar': './additional_files/include_regions_SRF.dat',
        'exptime_coadd': 2600,
        'rebinning_step': 0.010
    },
    'instruments': {
        'wavelength_step': 0.01000000,
        'use_rv_from_ccf': False,
        'use_analytical_rvs': False,
        'linear_fit_method': 'linear_curve_fit',
        'orders': None,
        'telluric': None,
        'telluric_template': None
    },
    'refraction': {
        'approach': 'individual_order',
        'method': 'polynomial',
        'fit_order': 5,
        'fit_iters': 3,
        'fit_sigma': 5,
        'knots_spacing': 25.00,
        'reference_night': False,
        'reference_instrument': False
    },
    'master-out': {
        'wavelength_range': [5880.000, 5906.000],
        'wavelength_step': 0.010000000,
        'binning_factor': 20,
        'use_smoothed': False,
        'use_composite': False,
        'boxcar_smoothing': 3,
    },
    'pca_parameters': {
        'iterations': 5,
        'ref_iteration': 0,
    },
    'settings': {
        'full_output': False
    }
}

""" List ok keywords that are copied from the instrument dictionary
    when they are not explicitily specified
"""
copy_from_instrument = [
    'telluric_template',
    'telluric',
    'mask',
    'use_rv_from_ccf',
    'use_analytical_rvs',
    'telescope',
    'spectral_selection',
    'apply_ESO_telluric_correction',
    'use_ESO_sky_correction',
    'use_ESO_deblazed',
    'savgol_window',
    'savgol_polyorder',
    'savgol_window_kms'
]
