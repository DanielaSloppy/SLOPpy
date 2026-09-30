from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *

from scipy.interpolate import UnivariateSpline

__all__ = ['compute_quick_transmission']

def compute_quick_transmission(config_in, lines_label):

    subroutine_name = 'quick_transmission'

    night_dict = from_config_get_nights(config_in, phase='transit')
    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]
    shared_data = load_from_cpickle('shared', config_in['output'])

    shared_selection = (shared_data['coadd']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['coadd']['wave'] < lines_dict['range'][1])
    binned_selection = (shared_data['binned']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['binned']['wave'] < lines_dict['range'][1])

    transmission_shared= {
        'subroutine': subroutine_name,
        'binned_wave': shared_data['binned']['wave'][binned_selection],
        'binned_step': shared_data['binned']['step'][binned_selection],
        'binned_size': int(np.sum(binned_selection))
    }

    import matplotlib.pyplot as plt

    for night in night_dict:

        #try:
        #    preparation = load_from_cpickle('transmission_preparation',
        #                                  config_in['output'],
        #                                  night)
        #    print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
        #    continue
        #except:
        #    print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
        #    print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        if config_in['master-out'].get('use_composite', False):
            master_out = load_from_cpickle('master_out_composite', config_in['output'], night)
            print('  Using composite master-out from all nights')
        else:
            master_out = load_from_cpickle('master_out', config_in['output'], night)

        if config_in['master-out'].get('use_smoothed', False):
            master_out['rescaled'] = master_out['smoothed']
            master_out['rescaled_err'] = master_out['smoothed_err']
            print('  Using smoothed master-out')

        quick_transmission = {
            'subroutine': subroutine_name,
        }

        first_obs = lists['observations'][0]

        quick_transmission['n_pixels'] = input_data[first_obs]['n_pixels']
        quick_transmission['n_orders'] = input_data[first_obs]['n_orders']
        quick_transmission['wave'] = input_data[first_obs]['wave']
        quick_transmission['step'] = input_data[first_obs]['step']
        master_raw = np.zeros([quick_transmission['n_orders'], quick_transmission['n_pixels']])


        for obs in lists['transit_out']:
            master_raw += input_data[obs]['e2ds']
        quick_transmission['master_raw'] = master_raw.copy()
        quick_transmission['master'] = master_raw / np.nanmedian(master_raw)

        blaze = np.ones([quick_transmission['n_orders'], quick_transmission['n_pixels']])


        for obs in lists['observations']:

            quick_transmission[obs] = {}

            transmission_raw = input_data[obs]['e2ds'] / quick_transmission['master']
            quick_transmission[obs]['transmission_raw'] = transmission_raw.copy()
            quick_transmission[obs]['transmission'] = transmission_raw / np.nanmedian(transmission_raw)
            quick_transmission[obs]['wave'] = input_data[obs]['wave'] #Added for plotting purpose only


            quick_transmission[obs]['binned'] = \
                rebin_2d_to_1d(quick_transmission['wave'],
                            quick_transmission['step'],
                            quick_transmission[obs]['transmission'],
                            blaze,
                            transmission_shared['binned_wave'],
                            transmission_shared['binned_step'],
                            rv_shift=0,
                            preserve_flux=False)

            #plt.scatter(quick_transmission['wave'],
            #            quick_transmission[obs]['transmission'],
            #            c='C1', s=1, zorder=3, alpha=0.25)

        for obs in lists['transit_out']:

            plt.scatter(transmission_shared['binned_wave'],
                        quick_transmission[obs]['binned'],
                        c='b', s=1, zorder=10, alpha=0.5)

        for obs in lists['transit_full']:

            plt.scatter(transmission_shared['binned_wave'],
                        quick_transmission[obs]['binned']-0.04,
                        c='r', s=1, zorder=10, alpha=0.5)

        plt.xlim(transmission_shared['binned_wave'][0], transmission_shared['binned_wave'][-1])
        plt.show()


    print()
    """ Keep going from here after preparation, unless the subroutines has been called just
        to preform the data preparation step
    """

