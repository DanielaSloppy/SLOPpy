from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from astropy.convolution import convolve, Box1DKernel

__all__ = ['compute_master_out', 'plot_master_out', 'plot_compare_master_out']
subroutine_name = 'master_out'


def compute_master_out(config_in):

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    master_out_composite = {
        'subroutine': 'master_out',
        'wave': shared_data['coadd']['wave'],
        'step': shared_data['coadd']['step'],
        'size': shared_data['coadd']['size'],
    }
    """ The composite master-out is computed separately for transit and eclipse nights,
        so that the master-out of eclipse nights only includes spectra with the planet
        behind the star, and vice versa """
    wmean_wflux = {phase: np.zeros(master_out_composite['size']) for phase in ['transit', 'eclipse']}
    wmean_weight = {phase: np.zeros(master_out_composite['size']) for phase in ['transit', 'eclipse']}

    box_kernel = Box1DKernel(config_in['master-out'].get('boxcar_smoothing', 1))

    for night in night_dict:

        phase = get_night_phase(night_dict, night)

        try:
            master_out = load_from_cpickle('master_out', config_in['output'], night)
            wmean_wflux[phase] += master_out['SRF']['rescaled'] / master_out['SRF']['rescaled_err'] ** 2
            wmean_weight[phase] += 1. / master_out['SRF']['rescaled_err'] ** 2
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': 'master_out'
        }

        master_out = {
            'subroutine': subroutine_name,
            'wave': shared_data['coadd']['wave'],
            'step': shared_data['coadd']['step'],
            'size': shared_data['coadd']['size'],
            'total_flux': np.zeros(shared_data['coadd']['size'], dtype=np.double),
            'total_flux_err': np.zeros(shared_data['coadd']['size'], dtype=np.double)
        }

        for obs in lists['transit_out']:

            processed[obs] = {}

            processed[obs]['rescaling'], \
            processed[obs]['rescaled'], \
            processed[obs]['rebinned_err'] = perform_rescaling(
                input_data[obs]['wave'], input_data[obs]['e2ds'], input_data[obs]['e2ds_err'],
                observational_pams['wavelength_rescaling'])

            preserve_flux = input_data[obs].get('absolute_flux', True)

            processed[obs]['rebinned'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds'],
                               calib_data['blaze'],
                               master_out['wave'],
                               master_out['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'])

            processed[obs]['rebinned_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds_err'],
                               calib_data['blaze'],
                               master_out['wave'],
                               master_out['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'])

            master_out['total_flux'] += processed[obs]['rebinned']
            master_out['total_flux_err'] += processed[obs]['rebinned_err']**2.0

        master_out['total_flux_err'] = np.sqrt(master_out['total_flux_err'])

        master_out['rescaling'], \
        master_out['rescaled'], \
        master_out['rescaled_err'] = perform_rescaling(
            master_out['wave'], master_out['total_flux'],
            master_out['total_flux_err'],
            observational_pams['wavelength_rescaling'])

        master_out['rescaled'], master_out['rescaled_err'], master_out['null'] = \
            replace_values_errors(master_out['rescaled'], master_out['rescaled_err'],
                                  threshold=0.0001, replacement=1.0000)

        master_out['smoothed'] = convolve(master_out['rescaled'].copy(), box_kernel)
        master_out['smoothed_err'] = np.sqrt(convolve((master_out['rescaled_err'])**2, box_kernel))

        sel = (master_out['smoothed']<0.01) | (master_out['smoothed']>1.5)
        master_out['smoothed'][sel] = 1.0
        master_out['smoothed_err'][sel] = 1.0

        selection = (master_out['wave']>0)
        spline_iter = 5
        for n_iter in range(0, spline_iter):

            residuals = master_out['rescaled'] / master_out['smoothed']
            wave = master_out_composite['wave']

            """ picking the number of knots """
            nknots = ((np.amax(wave) - np.amin(wave)) / config_in['master-out'].get('spline_step', 0.10))
            """ picking the indices of the knots"""
            idx_knots = (np.arange(1, len(wave[selection]) - 1, (len(wave[selection]) - 2.) / nknots)).astype('int')
            """ passing from indices to knots values """
            knots = wave[selection][idx_knots]

            coeff = sci_int.splrep(wave[selection], residuals[selection], task=-1, k=2, t=knots)
            spline = sci_int.splev(wave, coeff)

            dif = residuals - spline
            std = np.std(dif)
            selection = np.where(np.abs(dif) < 4 * std)  # & (refraction[obs]['flag'])

        master_out['spline'] = spline
        master_out['smoothed'] *= spline
        master_out['smoothed_err'] *= spline

        master_out['SRF'] = {}

        master_out['SRF']['rescaled']= \
            rebin_1d_to_1d(master_out['wave'],
                           master_out['step'],
                           master_out['rescaled'],
                           master_out['wave'], master_out['step'],
                           rv_shift=observational_pams['rv_shift_ORF2SRF_res'],
                           preserve_flux=False)

        master_out['SRF']['rescaled_err']= \
            rebin_1d_to_1d(master_out['wave'],
                           master_out['step'],
                           master_out['rescaled_err'],
                           master_out['wave'], master_out['step'],
                           rv_shift=observational_pams['rv_shift_ORF2SRF_res'],
                           preserve_flux=False,
                           is_error=True)

        wmean_wflux[phase] += master_out['SRF']['rescaled']/master_out['SRF']['rescaled_err']**2
        wmean_weight[phase] += 1./master_out['SRF']['rescaled_err']**2



        """
        rv_shift = observational_pams['BERV_avg'] - observational_pams['RV_star']['intercept']
        # bringing the master-out to the aboslute reference system 
        wave_shifted, _ = shift_wavelength(master_out['wave'], master_out['step'], rv_shift)

        # master-out is printed to .dat for compatibility with other programs
        master_data_out = get_filename('master_out', config_in['output'], night, extension=".dat")
        file_out = open(master_data_out, 'w')
        for w, f, e in  zip(wave_shifted, master_out['rescaled'], master_out['rescaled_err']):
            file_out.write('{0:10.4f} {1:f} {2:f}\n'.format(w,f,e))
        file_out.close()
        print()
        print("NON OPTIMAL MASTER-OUT DAT FILE!!!!")
        print()
        """

        save_to_cpickle('master_out_processed', processed, config_in['output'], night)
        save_to_cpickle('master_out', master_out, config_in['output'], night)

    composite_SRF = {}
    for phase in ['transit', 'eclipse']:
        """ skip the phase if no night of that type has been analyzed """
        if not np.any(wmean_weight[phase] > 0.):
            continue
        with np.errstate(divide='ignore', invalid='ignore'):
            composite_SRF[phase] = {
                'rescaled': wmean_wflux[phase]/wmean_weight[phase],
                'rescaled_err': np.sqrt(1./wmean_weight[phase])
            }
        composite_SRF[phase]['smoothed'] = convolve(composite_SRF[phase]['rescaled'].copy(), box_kernel)
        composite_SRF[phase]['smoothed_err'] = \
            np.sqrt(convolve((composite_SRF[phase]['rescaled_err']) ** 2, box_kernel))

    print()
    for night in night_dict:
        phase = get_night_phase(night_dict, night)
        try:
            master_out_composite = load_from_cpickle('master_out_composite', config_in['output'], night)
            """ recompute if the stored composite master-out has been obtained for the other phase;
                composites saved without the phase keyword were obtained from transit nights """
            if master_out_composite.get('phase', 'transit') != phase:
                raise ValueError('Stored composite master-out computed with different settings')
            print("{0:45s} Night:{1:15s}   {2:s}".format('master_out_composite', night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format('master_out_composite', night, 'Computing'))
            print()

        master_out_composite['SRF'] = composite_SRF[phase]
        master_out_composite['phase'] = phase

        master_out = load_from_cpickle('master_out', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        master_out_composite['rescaled']= \
            rebin_1d_to_1d(master_out_composite['wave'],
                           master_out_composite['step'],
                           master_out_composite['SRF']['rescaled'],
                           master_out_composite['wave'],
                           master_out_composite['step'],
                           rv_shift=-observational_pams['rv_shift_ORF2SRF_res'],
                           preserve_flux=False)

        master_out_composite['rescaled_err']= \
            rebin_1d_to_1d(master_out_composite['wave'],
                           master_out_composite['step'],
                           master_out_composite['SRF']['rescaled_err'],
                           master_out_composite['wave'],
                           master_out_composite['step'],
                           rv_shift=-observational_pams['rv_shift_ORF2SRF_res'],
                           preserve_flux=False,
                           is_error=True)

        master_out_composite['smoothed'] = \
            rebin_1d_to_1d(master_out_composite['wave'],
                           master_out_composite['step'],
                           master_out_composite['SRF']['smoothed'],
                           master_out_composite['wave'],
                           master_out_composite['step'],
                           rv_shift=-observational_pams['rv_shift_ORF2SRF_res'],
                           preserve_flux=False)

        master_out_composite['smoothed_err']= \
            rebin_1d_to_1d(master_out_composite['wave'],
                           master_out_composite['step'],
                           master_out_composite['SRF']['smoothed_err'],
                           master_out_composite['wave'],
                           master_out_composite['step'],
                           rv_shift=-observational_pams['rv_shift_ORF2SRF_res'],
                           preserve_flux=False,
                           is_error=True)

        #master_out_composite['smoothed'] = convolve(master_out_composite['rescaled'].copy(), box_kernel)
        #master_out_composite['smoothed_err'] = \
        #    np.sqrt(convolve((master_out_composite['rescaled_err']) ** 2, box_kernel))

        sel = (master_out_composite['smoothed']<0.01) | (master_out_composite['smoothed']>1.5)
        master_out_composite['smoothed'][sel] = 1.0
        master_out_composite['smoothed_err'][sel] = 1.0

        selection = (master_out_composite['wave']>0)
        spline_iter = 5
        for n_iter in range(0, spline_iter):

            residuals = master_out['rescaled'] / master_out_composite['smoothed']
            wave = master_out_composite['wave']

            """ picking the number of knots """
            nknots = ((np.amax(wave) - np.amin(wave)) / config_in['master-out'].get('spline_step', 0.10))
            """ picking the indices of the knots"""
            idx_knots = (np.arange(1, len(wave[selection]) - 1, (len(wave[selection]) - 2.) / nknots)).astype('int')
            """ passing from indices to knots values """
            knots = wave[selection][idx_knots]

            coeff = sci_int.splrep(wave[selection], residuals[selection], task=-1, k=2, t=knots)
            spline = sci_int.splev(wave, coeff)

            dif = residuals - spline
            std = np.std(dif)
            selection = np.where(np.abs(dif) < 4 * std)  # & (refraction[obs]['flag'])

        master_out_composite['spline'] = spline
        master_out_composite['smoothed'] *= spline
        master_out_composite['smoothed_err'] *= spline

        save_to_cpickle('master_out_composite', master_out_composite, config_in['output'], night)


def plot_master_out(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)

    if night_input=='':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_master_out                            Night: ", night)

        """ Retrieving the analysis"""
        try:
            master_out = load_from_cpickle('master_out', config_in['output'], night)
            master_out_composite = load_from_cpickle('master_out_composite', config_in['output'], night)
        except:
            print()
            print("No master_out , no plots")
            continue

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        plt.figure(figsize=(12, 6))
        plt.title('Master out - night ' + night)

        lists = load_from_cpickle('lists', config_in['output'], night)
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)


        obs = lists['transit_out'][0]
        plt.scatter(input_data[obs]['wave'],
                 calib_data['blaze'],
                 color='C3', zorder=3., label='blaze', alpha=0.25)

        plt.errorbar(master_out['wave'],
                     master_out['rescaled'],
                     yerr=master_out['rescaled_err'],
                     fmt='.', c='C0', label='master-out ' + night)
        plt.plot(master_out['wave'],
                 master_out['smoothed'],
                 color='C1', zorder=3., label='smoothed master-out ' + night)
        plt.scatter(master_out_composite['wave'],
                    master_out_composite['rescaled'],
                    s=2, c='C3')
        plt.plot(master_out_composite['wave'],
                 master_out_composite['smoothed'],
                 c='C3', label='composite master-out')
        plt.scatter(master_out['wave'],
                    master_out['rescaled']/master_out['smoothed']*master_out['spline']+0.05,
                    s=2, c='C4', label='rescaled/smoothed')
        plt.scatter(master_out['wave'],
                    master_out['rescaled']/master_out_composite['smoothed']*master_out_composite['spline']+0.1,
                    s=2, c='C5', label='rescaled/ comp smoothed')
        plt.plot(master_out['wave'],
                 master_out['spline']+0.05,
                 c='C7', label='spline fit of the residuals')
        plt.plot(master_out['wave'],
                 master_out_composite['spline']+0.1,
                 c='C7')
        plt.ylim(0, 1.25)
        plt.xlabel('$\lambda$ [$\AA$]')
        plt.ylabel('Rescaled flux')
        plt.legend()
        plt.show()


def plot_compare_master_out(config_in):

    plt.figure(figsize=(12, 6))
    plt.title('Master out - comparison between nights ')
    night_dict = from_config_get_nights(config_in)

    for i, night in enumerate(night_dict):

        """ Retrieving the analysis"""
        try:
            master_out = load_from_cpickle('master_out', config_in['output'], night)
            master_out_composite = load_from_cpickle('master_out_composite', config_in['output'], night)
        except:
            continue

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        plt.errorbar(master_out['wave'],
                     master_out['SRF']['rescaled'],
                     yerr=master_out['SRF']['rescaled_err'],
                     fmt='.', c='C'+repr(i), label='master-out ' + night, alpha=0.5)

        if i == 0:
            plt.plot(master_out_composite['wave'],
                     master_out_composite['SRF']['rescaled'],
                     color='k', zorder=10, label='composite master-out')

    plt.ylim(0, 1.25)
    plt.xlabel('$\lambda$ [$\AA$]')
    plt.ylabel('Rescaled flux')
    plt.legend()
    plt.show()