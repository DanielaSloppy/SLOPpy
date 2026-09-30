from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.shortcuts import *

""" Subroutines shared by transmission_spectrum_preparation and emission_spectrum_preparation """

__all__ = ['load_master_out_for_preparation',
           'rebin_master_out_to_observation']


def load_master_out_for_preparation(config_in, night):
    """
    Load the master-out of the night, following the options of the master-out section
    :param config_in: configuration dictionary
    :param night: name of the night
    :return: master-out dictionary
    """
    if config_in['master-out'].get('use_composite', False):
        master_out = load_from_cpickle('master_out_composite', config_in['output'], night)
        print('  Using composite master-out from all nights')
    else:
        master_out = load_from_cpickle('master_out', config_in['output'], night)

    if config_in['master-out'].get('use_smoothed', False):
        master_out['rescaled'] = master_out['smoothed']
        master_out['rescaled_err'] = master_out['smoothed_err']
        print('  Using smoothed master-out')

    return master_out


def rebin_master_out_to_observation(master_out, input_data_obs, rv_shift_ORF2SRF, n_orders):
    """
    Bring back the 1D master-out to the ORF and rebin it to the 2D scale of the observation
    :param master_out: master-out dictionary
    :param input_data_obs: dictionary of the observation
    :param rv_shift_ORF2SRF: RV shift from the observer to the stellar rest frame
    :param n_orders: number of orders
    :return: rebinned master-out and its error
    """
    rebinned = rebin_1d_to_2d(master_out['wave'],
                              master_out['step'],
                              master_out['rescaled'],
                              input_data_obs['wave'],
                              input_data_obs['step'],
                              rv_shift=-rv_shift_ORF2SRF,
                              preserve_flux=False)

    rebinned_err = rebin_1d_to_2d(master_out['wave'],
                                  master_out['step'],
                                  master_out['rescaled_err'],
                                  input_data_obs['wave'],
                                  input_data_obs['step'],
                                  rv_shift=-rv_shift_ORF2SRF,
                                  preserve_flux=False,
                                  is_error=True)

    for order in range(0, n_orders):
        rebinned[order, :], rebinned_err[order, :], _ = \
            replace_values_errors_with_interpolation_1d(rebinned[order, :],
                                                        rebinned_err[order, :],
                                                        less_than=0.001, greater_than=5.0000)

    return rebinned, rebinned_err
