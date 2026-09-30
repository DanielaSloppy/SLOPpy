from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.pca_preparation import compute_pca_preparation

from PyAstronomy import pyasl

from scipy.interpolate import UnivariateSpline


__all__ = ['compute_sysrem_correction',
           'plot_sysrem_correction']

def compute_sysrem_correction(config_in):

    subroutine_name = 'sysrem_correction'
    compute_pca_preparation(config_in)
    print()

    night_dict = from_config_get_nights(config_in, phase='transit')
    pca_parameters = from_config_get_pca_parameters(config_in)

    for night in night_dict:

        try:
            sysrem_output = load_from_cpickle('transmission_preparation',
                                          config_in['output'],
                                          night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()


        n_iter = pca_parameters.get('iterations',5 )
        ref_iter = pca_parameters.get('ref_iteration',3 )
        sysrem_output = {
            'subroutine': subroutine_name,
            'iterations': n_iter,
            'pca_output': True,
            'ref_iteration': ref_iter
        }


        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        preparation = load_from_cpickle('pca_preparation', config_in['output'], night)

        n_obs, n_orders, n_pixels = np.shape(preparation['stack_e2ds']) 

        model_iter = np.ones([n_iter, n_obs, n_orders, n_pixels], dtype=np.double)
        model_out = np.ones([n_iter, n_obs, n_orders, n_pixels], dtype=np.double)

        for order in range(0, observational_pams['n_orders']):
            obs = preparation['stack_e2ds'][:,order,:]
            sigs = preparation['stack_e2ds_err'][:,order,:]

            sr = pyasl.SysRem(obs, sigs)
            previous_residuals = obs.copy()

            for it in range(0, n_iter):
                r, a, c = sr.iterate()
                model_iter[it,:,order,:] =  previous_residuals - r
                previous_residuals = r.copy()


        for it in range(0, n_iter):

            # Full model is the sum of all the models until the given iteration 
            model_out[it,:, :, :] = np.sum(model_iter[:it+1,:, :, :], axis=0)

            #import matplotlib.pyplot as plt
            #plt.figure()
            #plt.title("Model " + repr(it))
            #plt.imshow( model_out[it,:, 10, :], origin='lower', aspect="auto")
            #plt.show()

            it_string = str(it).zfill(2)

            sysrem_output[it_string] = {
                'model': model_out[it,:, :, :]
            }

            for i_obs, obs in enumerate(lists['observations']):

                sysrem_output[it_string][obs] = {}
                sysrem_output[it_string][obs]['ratio']  = preparation['stack_e2ds'][i_obs, :, :]/model_out[it, i_obs, :, :]
                sysrem_output[it_string][obs]['ratio_err']  = preparation['stack_e2ds_err'][i_obs, :, :]/model_out[it, i_obs, :, :]


        save_to_cpickle('transmission_preparation', sysrem_output, config_in['output'], night)

    print()
    """ Keep going from here after preparation, unless the subroutines has been called just
        to preform the data preparation step
    """


def plot_sysrem_correction(config_in, night_input=''):

    subroutine_name = 'transmission_spectrum_preparation'

    night_dict = from_config_get_nights(config_in, phase='transit')

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)




    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        # ! To be removed when testing is done
        # ! This plots do not make any sense anymore 
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])


        """ Retrieving the analysis"""
        try:
            preparation = load_from_cpickle('transmission_preparation', config_in['output'], night)
        except:
            print("No transmission spectrum results, no plots")
            print()
            continue


        #from SLOPpy.subroutines.lines_fit_functions import logprob_case12
        from matplotlib.colors import BoundaryNorm
        from matplotlib.ticker import MaxNLocator

        len_y = len(lists['observations'])
        len_x = 4096
        order= 11

        time_from_transit = np.empty(len_y, dtype=np.double)
        plot_data = np.empty([len_y, len_x], dtype=np.double)

        for i_obs, obs in enumerate(lists['observations']):
            print(np.median(preparation[obs]['deblazed'][order ,:]))
            time_from_transit[i_obs] =  input_data[obs]['BJD'] - observational_pams['time_of_transit']
            plot_data[i_obs, :] = preparation[obs]['deblazed'][order ,:]/ np.median(preparation[obs]['ratio'][order ,:])
            wave = preparation[obs]['wave'][order, :]

        wave_meshgrid, time_meshgrid = np.meshgrid(wave, time_from_transit)

        print('COOLWARM')
        cmap = plt.get_cmap('coolwarm')

        levels = MaxNLocator(nbins=15).tick_values(
            plot_data.min(), plot_data.max())
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plt.figure(figsize=(15, 10))

        PCF = plt.contourf(wave_meshgrid, time_meshgrid,
                            plot_data, levels=levels, cmap=cmap)
        cbar = plt.colorbar(PCF)
        cbar.ax.set_ylabel('Intensity')
        plt.show()

        res = plot_data * 1.
        from scipy.interpolate import UnivariateSpline
        for ii in range(0,4096):
            spl = UnivariateSpline(time_from_transit, plot_data[:, ii])
            val = spl(time_from_transit)
            res[:,ii] -= val
            res[:,ii] /= val



        print('COOLWARM')
        cmap = plt.get_cmap('coolwarm')

        levels = MaxNLocator(nbins=10).tick_values(
            -0.05, 0.05)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plt.figure(figsize=(15, 10))

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
        #ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for obs in lists['transit_in']:

            preparation[obs]['rescaling'], \
            preparation[obs]['rescaled'], \
            preparation[obs]['rescaled_err'] = perform_rescaling(
                preparation[obs]['wave'],
                preparation[obs]['deblazed'] / (input_data[obs]['step'] / np.median(input_data[obs]['step'])),
                preparation[obs]['deblazed_err'] / (input_data[obs]['step'] / np.median(input_data[obs]['step'])),
                observational_pams['wavelength_rescaling'])

            ax1.scatter(preparation[obs]['wave'],
                    preparation[obs]['rescaled'],
                    s=1, alpha=0.25,
                    color=colors_plot['mBJD'][obs])

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
