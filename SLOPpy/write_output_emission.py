from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.emission_spectrum import compute_emission_spectrum_core, plot_emission_spectrum_core

__all__ = ['write_output_emission', 'plot_output_emission',
           'write_output_emission_planetRF', 'plot_output_emission_planetRF',
           'write_output_emission_stellarRF', 'plot_output_emission_stellarRF',
           'write_output_emission_observerRF', 'plot_output_emission_observerRF']

subroutine_name = 'write_output_emission'


def write_output_emission_planetRF(config_in):
    write_output_emission(config_in, reference='planetRF')


def plot_output_emission_planetRF(config_in, night_input=''):
    plot_output_emission(config_in, night_input, reference='planetRF')


def write_output_emission_stellarRF(config_in):
    write_output_emission(config_in, reference='stellarRF')


def plot_output_emission_stellarRF(config_in, night_input=''):
    plot_output_emission(config_in, night_input, reference='stellarRF')


def write_output_emission_observerRF(config_in):
    write_output_emission(config_in, reference='observerRF')


def plot_output_emission_observerRF(config_in, night_input=''):
    plot_output_emission(config_in, night_input, reference='observerRF')


def write_output_emission(config_in, reference='planetRF'):
    """
    Emission spectrum of each eclipse night over the full spectrum, or over the range
    specified in the full_spectrum section of the configuration file.
    Unlike write_output_transmission, no CLV+RM correction and no continuum normalization are
    applied, as the emission spectra have already been normalized by emission_spectrum_preparation
    """
    fullspectrum_dict = from_config_get_fullspectrum_parameters(config_in)

    compute_emission_spectrum_core(config_in, subroutine_name, fullspectrum_dict.get('range', None),
                                   reference=reference)


def plot_output_emission(config_in, night_input='', reference='planetRF'):

    fullspectrum_dict = from_config_get_fullspectrum_parameters(config_in)
    shared_data = load_from_cpickle('shared', config_in['output'])

    plot_range = fullspectrum_dict.get('plot_range',
                                       fullspectrum_dict.get('range', shared_data['coadd']['wavelength_range']))

    plot_emission_spectrum_core(config_in, subroutine_name, plot_range,
                                reference=reference, night_input=night_input)
