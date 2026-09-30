from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.constants import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.preparation_subroutines import *

from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter


__all__ = ['compute_emission_spectrum_preparation',
           'plot_emission_spectrum_preparation',
           'compute_savgol_continuum']


def compute_savgol_continuum(flux, wave, step, window=301, polyorder=2, window_kms=None):
    """
    Continuum of each order of a 2D spectrum, estimated with a Savitzky-Golay filter
    :param flux: 2D array (n_orders x n_pixels)
    :param wave: 2D wavelength array
    :param step: 2D wavelength step array
    :param window: width of the filter window in pixels
    :param polyorder: order of the polynomial fitted within the window
    :param window_kms: width of the filter window in km/s. When provided, it overrides
        window and it is converted into pixels order by order, using the local
        velocity step, so that the same value can be used for any instrument
    :return: 2D array with the continuum, with zero values replaced by unity
    """
    n_orders, n_pixels = np.shape(flux)
    continuum = np.empty([n_orders, n_pixels], dtype=np.double)

    """ the window must be odd, larger than polyorder and not larger than the order """
    max_window = n_pixels if n_pixels % 2 == 1 else n_pixels - 1
    min_window = polyorder + 2 if polyorder % 2 == 1 else polyorder + 1

    for order in range(0, n_orders):
        if window_kms is not None:
            velocity_step = speed_of_light_km * np.median(step[order, :] / wave[order, :])
            window_order = int(np.round(window_kms / velocity_step))
        else:
            window_order = int(window)

        if window_order % 2 == 0:
            window_order += 1
        window_order = min(max(window_order, min_window), max_window)

        continuum[order, :] = savgol_filter(flux[order, :], window_order, polyorder)

    continuum[continuum == 0] = 1
    return continuum


def compute_emission_spectrum_preparation(config_in):

    subroutine_name = 'emission_spectrum_preparation'

    """ only nights with phase: eclipse are analyzed here """
    night_dict = from_config_get_nights(config_in, phase='eclipse')

    for night in night_dict:

        savgol_pams = {
            'window': int(night_dict[night].get('savgol_window', 301)),
            'polyorder': int(night_dict[night].get('savgol_polyorder', 2)),
            'window_kms': night_dict[night].get('savgol_window_kms', None)
        }

        try:
            preparation = load_from_cpickle('emission_preparation',
                                          config_in['output'],
                                          night)
            """ recompute if the stored spectra were obtained with a different normalization """
            if preparation.get('savgol_pams', None) != savgol_pams:
                raise ValueError('Stored preparation computed with different settings')
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        master_out = load_master_out_for_preparation(config_in, night)

        if savgol_pams['window_kms'] is not None:
            print('  Continuum normalization: Savitzky-Golay filter, window {0:.1f} km/s, polyorder {1:d}'.format(
                savgol_pams['window_kms'], savgol_pams['polyorder']))
        else:
            print('  Continuum normalization: Savitzky-Golay filter, window {0:d} pixels, polyorder {1:d}'.format(
                savgol_pams['window'], savgol_pams['polyorder']))

        preparation = {
            'subroutine': subroutine_name,
            'savgol_pams': savgol_pams,
        }

        for obs in lists['observations']:

            preparation[obs] = {}

            preparation[obs]['master_out'] = {}
            preparation[obs]['wave'] = input_data[obs]['wave'] #Added for plotting purpose only

            """ Step 1+2): bring back the master-out to the ORF and rebin the 1D master-out to the 2D observation scale"""
            preparation[obs]['master_out']['rebinned'], \
            preparation[obs]['master_out']['rebinned_err'] = \
                rebin_master_out_to_observation(master_out,
                                                input_data[obs],
                                                observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                                                observational_pams['n_orders'])

            """ Step 3): obtain the emission spectrum for this observation, by subtracting the
                master-out from the observation after both have been normalized by their continuum.
                The result is stored under the same keywords used for the transmission spectrum
            """
            e2ds_cont = compute_savgol_continuum(input_data[obs]['e2ds'],
                                                 input_data[obs]['wave'],
                                                 input_data[obs]['step'],
                                                 **savgol_pams)
            master_out_cont = compute_savgol_continuum(preparation[obs]['master_out']['rebinned'],
                                                       input_data[obs]['wave'],
                                                       input_data[obs]['step'],
                                                       **savgol_pams)

            preparation[obs]['ratio'] = input_data[obs]['e2ds'] / e2ds_cont - \
                                      preparation[obs]['master_out']['rebinned'] / master_out_cont
            preparation[obs]['ratio_err'] = np.sqrt((input_data[obs]['e2ds_err'] / e2ds_cont)**2 +
                                                  (preparation[obs]['master_out']['rebinned_err'] /
                                                   master_out_cont)**2)

            preparation[obs]['ratio_precleaning'] = preparation[obs]['ratio'].copy()
            preparation[obs]['ratio_precleaning_err'] = preparation[obs]['ratio_err'].copy()


        if night_dict[night].get('spline_residuals', True):

            print()
            print('   Cleaning for telluric residuals with Univariate Spline - threshold about 0.05')
            """ emission spectra are centred on zero: they are shifted to unity before the
                spline fit, instead of being normalized by their median """
            for order in range(0, observational_pams['n_orders']):
                obs_reference =  lists['observations'][0]

                len_y = len(lists['observations'])
                len_x = len(preparation[obs_reference]['wave'][order, :])

                time_from_transit = np.empty(len_y, dtype=np.double)
                data_array = np.empty([len_y, len_x], dtype=np.double)
                for i_obs, obs in enumerate(lists['observations']):
                    time_from_transit[i_obs] =  input_data[obs]['BJD'] - observational_pams['time_of_transit']
                    data_array[i_obs, :] = preparation[obs]['ratio_precleaning'][order ,:] + 1.

                res = data_array * 1.
                val = np.empty([len_y, len_x], dtype=np.double)

                for ii in range(0, len_x):
                    spl = UnivariateSpline(time_from_transit, data_array[:, ii])
                    val[:,ii] = spl(time_from_transit)
                    res[:,ii] -= val[:,ii]
                    res[:,ii] /= val[:,ii]

                sel = np.abs(res) > 0.05

                for i_obs, obs in enumerate(lists['observations']):
                    if np.sum(sel[i_obs]) > 0:
                        preparation[obs]['ratio'][order, sel[i_obs]] = val[i_obs, sel[i_obs]] - 1.
                        preparation[obs]['ratio_err'][order, sel[i_obs]] *= 10.
        else:
            print()
            print('   Cleaning for telluric residuals NOT performed')


        for obs in lists['observations']:

            """ blaze and pixel-step dependence have already been removed by the continuum normalization """
            preparation[obs]['deblazed'] = preparation[obs]['ratio'].copy()
            preparation[obs]['deblazed_err'] = preparation[obs]['ratio_err'].copy()

            if not config_in['settings'].get('full_output', False):
                del preparation[obs]['master_out']

        save_to_cpickle('emission_preparation', preparation, config_in['output'], night)

    print()


def plot_emission_spectrum_preparation(config_in, night_input=''):

    subroutine_name = 'emission_spectrum_preparation'

    night_dict = from_config_get_nights(config_in, phase='eclipse')

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])

        """ Retrieving the analysis"""
        try:
            preparation = load_from_cpickle('emission_preparation', config_in['output'], night)
        except:
            print("No emission spectrum results, no plots")
            print()
            continue

        from matplotlib.colors import BoundaryNorm
        from matplotlib.ticker import MaxNLocator

        """ emission spectra are centred on zero: they are shifted to unity for plotting purposes """
        obs_reference = lists['observations'][0]
        len_y = len(lists['observations'])
        len_x = np.shape(preparation[obs_reference]['deblazed'])[1]
        order = 11

        time_from_transit = np.empty(len_y, dtype=np.double)
        plot_data = np.empty([len_y, len_x], dtype=np.double)

        for i_obs, obs in enumerate(lists['observations']):
            time_from_transit[i_obs] =  input_data[obs]['BJD'] - observational_pams['time_of_transit']
            plot_data[i_obs, :] = preparation[obs]['deblazed'][order ,:] + 1.
            wave = preparation[obs]['wave'][order, :]

        wave_meshgrid, time_meshgrid = np.meshgrid(wave, time_from_transit)

        cmap = plt.get_cmap('coolwarm')

        levels = MaxNLocator(nbins=21).tick_values(0.90, 1.10)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        plt.figure(figsize=(15, 10))
        plt.title('Emission map in observer reference frame (shifted by +1)\n {0:s}'.format(night))

        PCF = plt.contourf(wave_meshgrid, time_meshgrid,
                            plot_data, levels=levels, cmap=cmap)
        cbar = plt.colorbar(PCF)
        cbar.ax.set_ylabel('Intensity')
        plt.show()

        if night_dict[night].get('spline_residuals', True):
            res = plot_data * 1.
            for ii in range(0, len_x):
                spl = UnivariateSpline(time_from_transit, plot_data[:, ii])
                val = spl(time_from_transit)
                res[:,ii] -= val
                res[:,ii] /= val

            cmap = plt.get_cmap('coolwarm')

            levels = MaxNLocator(nbins=10).tick_values(-0.05, 0.05)
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            plt.figure(figsize=(15, 10))
            plt.title('Residuals after dividing by UnivariateSpline Spline\n {0:s}'.format(night))

            PCF = plt.contourf(wave_meshgrid, time_meshgrid,
                                res, levels=levels, cmap=cmap)
            cbar = plt.colorbar(PCF)
            cbar.ax.set_ylabel('Intensity')
            plt.show()

        """ Creation of the color array, based on the BJD of the observations
        """
        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))

        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax1 = plt.subplot(gs[0, 0])

        ax1.set_ylim(-0.10, 0.10)
        cbax1 = plt.subplot(gs[:, 1])

        for obs in lists['transit_in']:
            ax1.scatter(preparation[obs]['wave'],
                    preparation[obs]['deblazed'],
                    s=1, alpha=0.25,
                    color=colors_plot['mBJD'][obs])

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
