from __future__ import print_function, division
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle

import os
from os import path
import oyaml as yaml
from SLOPpy.subroutines.object_parameters import StarParameters, PlanetParameters
from SLOPpy.config_default import *
from SLOPpy.subroutines.kepler_exo import kepler_eclipse_time_offset

__all__ = ["save_to_cpickle",
           "load_from_cpickle",
           "delete_cpickle",
           "check_existence_cpickle",
           "get_filename",
           "load_yaml_file",
           "pars_input",
           "yaml_parser",
           "get_filelists",
           "from_config_get_nights",
           "get_night_phase",
           "from_config_get_instrument",
           "from_config_get_system",
           "from_config_get_pipeline",
           "from_config_get_planet",
           "from_config_get_star",
           "from_config_get_clv_rm",
           "from_config_refraction",
           "from_config_get_interstellar_lines",
           "from_config_get_transmission_lightcurve",
           "from_config_get_transmission",
           "from_config_get_molecfit",
           "from_config_get_transmission_mcmc",
           "from_config_get_spectral_lines",
           "from_config_get_interactive_plots",
           "from_config_get_pca_parameters",
           "from_config_get_fullspectrum_parameters"]

accepted_extensions = ['.yaml', '.yml', '.conf', '.config', '.input', ]


def save_to_cpickle(fname, dictionary, output, night='', lines='', it_string='', temp_folder=''):

    output_file = get_filename(fname, output, night, lines, it_string, temp_folder=temp_folder)
    pickle.dump(dictionary, open(output_file, "wb"))


def load_from_cpickle(fname, output, night='', lines='', it_string='', temp_folder=''):

    output_file = get_filename(fname, output, night, lines, it_string, temp_folder=temp_folder)
    return pickle.load(open(output_file, "rb"))


def delete_cpickle(fname, output, night='', lines='', it_string='', temp_folder=''):

    output_file = get_filename(fname, output, night, lines, it_string, temp_folder=temp_folder)
    os.remove(output_file)

def check_existence_cpickle(fname, output, night='', lines='', it_string='', temp_folder=''):

    output_file = get_filename(fname, output, night, lines, it_string, temp_folder=temp_folder)
    return path.isfile(output_file)


def get_filename(fname, output, night, lines='', it_string='', extension=".p", temp_folder=''):

    if len(temp_folder) > 0:
        str_lines = temp_folder + '/'+ output
    else:
        str_lines = output

    for str_input in [lines, night, fname, it_string]:
        if len(str_input) > 0:
            str_lines += '_' + str_input
    return str_lines + extension

    if lines == '':
        if night == '':
            return output + '_' + fname + extension
        else:
            return output + '_' + night + "_" + fname + extension
    else:
        if night == '':
            return output + '_' + lines + '_' + fname + extension
        else:
            return output + '_' + lines + '_' + night + "_" + fname + extension


def load_yaml_file(file_conf):
    # shortcut for jupyter notebook plots
    config_in = yaml_parser(file_conf)
    return pars_input(config_in)

def yaml_parser(file_conf):
    stream = open(file_conf, 'r')

    try:
        config_in = yaml.load(stream, Loader=yaml.FullLoader)
    except AttributeError:
        config_in = yaml.load(stream)
        print(' Consider updating YAML')
    except:
        print(' Some error happened while reading the configuration file')
        quit()

    if 'output' not in config_in:

        for extension in accepted_extensions:
            if file_conf.find(extension) > 0:
                output_name = file_conf.replace(extension, "")
                continue

        config_in['output'] = output_name

    return config_in


def pars_input(config_in):
    config_in['system'] = {}

    if 'settings' not in config_in:
        config_in['settings'] = config_default['settings'].copy()
    else:
        for key, key_val in config_default['settings'].items():
            if key not in config_in['settings']:
                 config_in['settings'][key] = key_val


    for instrument in config_in['instruments']:

        for key, key_val in config_default['instruments'].items():
            if key not in config_in['instruments'][instrument]:
                 config_in['instruments'][instrument][key] = key_val

        """ create the refraction dictionary if not listed under the instrument section"""
        if 'refraction' not in config_in['instruments'][instrument]:
            config_in['instruments'][instrument]['refraction'] = {}

        """ when the refractions parameters are not explicitely specified in this section, they are either inherited
        from the top level dictionary or copied from the default dictionary """
        for key, key_val in config_default['refraction'].items():
            if key not in config_in['instruments'][instrument]['refraction']:
                try:
                    config_in['instruments'][instrument]['refraction'][key] = config_in['refraction'][key]
                except:
                    config_in['instruments'][instrument]['refraction'][key] = key_val

    if 'master-out' not in config_in:
        config_in['master-out'] = config_default['master-out'].copy()
    else:
        for key, key_val in config_default['master-out'].items():
            if key not in config_in['master-out']:
                 config_in['master-out'][key] = key_val

    if 'shared' not in config_in['instruments']:
        if 'wavelength_step' not in config_in['master-out']:
            config_in['instruments']['shared'] = {
                'wavelength_step': 0.0100
            }
        else:
            config_in['instruments']['shared'] = {
                'wavelength_step': config_in['master-out']['wavelength_step']
            }

    if 'molecfit' not in config_in:
        config_in['molecfit'] = config_default['molecfit'].copy()
    else:
        for key, key_val in config_default['molecfit'].items():
            if key not in config_in['molecfit']:
                 config_in['molecfit'][key] = key_val

    if 'molecfit' not in config_in:
        config_in['molecfit'] = config_default['molecfit'].copy()
    else:
        for key, key_val in config_default['molecfit'].items():
            if key not in config_in['molecfit']:
                 config_in['molecfit'][key] = key_val



    for night in config_in['nights']:

        instrument = config_in['nights'][night]['instrument']

        """ keywords are inherited from the instrument dictionary, when not explicitely specified"""
        for key in copy_from_instrument:
            if key not in config_in['nights'][night]:
                try:
                    config_in['nights'][night][key] = config_in['instruments'][instrument][key]
                except KeyError:
                    continue

        if 'refraction' not in config_in['nights'][night]:
            config_in['nights'][night]['refraction'] = config_in['instruments'][instrument]['refraction'].copy()
        else:
            for key, key_val in config_in['instruments'][instrument]['refraction'].items():
                if key not in config_in['nights'][night]['refraction']:
                    config_in['nights'][night]['refraction'][key] = key_val

        #if 'master_out_method' not in config_in['nights'][night]:
        #    config_in['nights'][night]['master_out_method'] = None

        if config_in['nights'][night]['use_analytical_rvs'] and 'RV_semiamplitude' not in config_in['star']:
            print(" Missing RV_semiamplitude keyword for the star, the value will be computed from the RVs ")
            config_in['nights'][night]['use_analytical_rvs'] = False

        """ for eclipse nights, time_of_transit can be replaced by time_of_eclipse: the transit
            preceding the eclipse is used for the computation of the planetary RVs """
        if get_night_phase(config_in['nights'], night) == 'eclipse' \
                and 'time_of_transit' not in config_in['nights'][night]:
            if 'time_of_eclipse' not in config_in['nights'][night]:
                raise ValueError(
                    'Night {0:s}: time_of_transit or time_of_eclipse must be specified'.format(str(night)))
            planet_dict = config_in['planet']
            if planet_dict.get('orbit', 'circular') == 'circular':
                eccentricity, omega_rad = 0.0, np.pi / 2.0
            else:
                eccentricity, omega_rad = planet_dict['eccentricity'][0], planet_dict['omega'][0] * np.pi / 180.0
            config_in['nights'][night]['time_of_transit'] = \
                np.atleast_1d(config_in['nights'][night]['time_of_eclipse'])[0] - \
                kepler_eclipse_time_offset(np.atleast_1d(planet_dict['period'])[0], eccentricity, omega_rad)

    """ OLD approach to compute the RV of the planet, left here because it may be useful in the future
    try:

        _dict_star = {'mass': None, 'radius': None, 'gamma': None}
        for key in config_in['star']:
            _dict_star[key] = config_in['star'][key]

        config_in['system']['star'] = StarParameters(
            mass=_dict_star['mass'],
            radius=_dict_star['radius'])

        config_in['system']['common'] = {'degree': False, 'n_planets': 0, 'planets_list': []}
        for key in config_in['planets']['common']:
            config_in['system']['common'][key] = config_in['planets']['common'][key]

        for key in config_in['planets']:

            if key not in ['common']:
                config_in['system'][key] = PlanetParameters()
                config_in['system'][key].put_reference_epoch(config_in['system']['common']['Tref'])
                config_in['system'][key].put_RVparameters(
                    P=config_in['planets'][key]['P'],
                    K=config_in['planets'][key]['K'],
                    f=config_in['planets'][key]['Tc'],
                    e=config_in['planets'][key]['e'],
                    o=config_in['planets'][key]['o'],
                    degree=config_in['system']['common']['degree'])
                config_in['system'][key].put_RVplanet(config_in['planets'][key]['K_planet'])
                config_in['system'][key].put_star(config_in['system']['star'])
                config_in['system']['common']['n_planets'] += 1
                config_in['system']['common']['planets_list'].extend([key])
    except:
        pass
    """

    return config_in


def get_filelists(night_selected):
    """

    :param night_selected: usually the night_dict[night] dictionary from the main program
    :return:
    """

    """ List files are supposed to be in the same directory of the yaml file,
        NOT on the archive directory: in this way it is possible to try different
        combinations of nights and files without making a mess in the archive """

    files_list = np.atleast_1d(np.genfromtxt(night_selected['all'], dtype=str))

    try:
        files_transit_out = np.atleast_1d(np.genfromtxt(night_selected['out_transit'], dtype=str))
        files_transit_in = np.atleast_1d(np.genfromtxt(night_selected['in_transit'], dtype=str))
        files_transit_full = np.atleast_1d(np.genfromtxt(night_selected['full_transit'], dtype=str))
    except (FileNotFoundError, IOError):
        files_transit_out = None
        files_transit_in = None
        files_transit_full = None

    try:
        files_telluric = np.atleast_1d(np.genfromtxt(night_selected['telluric_list'], dtype=str))
    except (FileNotFoundError, IOError):
        files_telluric = None


    if night_selected['telluric'] is not None:
        files_star_telluric = np.atleast_1d(np.genfromtxt(night_selected['star_telluric'], dtype=str))
    else:
        files_star_telluric = None

    return files_list, files_transit_out, files_transit_in, files_transit_full, files_telluric, files_star_telluric


def from_config_get_nights(config_in, phase=None):
    """
    This subroutine creates a shortcut to the night list
    :param config_in:
    :param phase: if 'transit' or 'eclipse', only the nights covering that phase are returned
    :return: dictionary
    """
    if phase is None:
        return config_in['nights']
    return {night: night_pams for night, night_pams in config_in['nights'].items()
            if get_night_phase(config_in['nights'], night) == phase}


def get_night_phase(night_dict, night):
    """
    This subroutine returns the orbital phase covered by the night, as specified
    by the 'phase' keyword in the night section of the configuration file:
    'transit' (transmission spectrum, default) or 'eclipse' (emission spectrum)
    :param night_dict: dictionary of the nights
    :param night: name of the night
    :return: string
    """
    phase = night_dict[night].get('phase', 'transit')
    if phase not in ['transit', 'eclipse']:
        raise ValueError("Night {0:s}: phase keyword must be 'transit' or 'eclipse', "
                         "'{1:s}' given".format(str(night), str(phase)))
    return phase


def from_config_get_instrument(config_in):
    """
    This subroutine creates a shortcut to the instrument list
    :param config_in:
    :return: dictionary
    """
    return config_in['instruments']


def from_config_refraction(config_in, night):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    return config_in['nights'][night]['refraction']

def from_config_get_transmission_lightcurve(config_in):
    """
    This subroutine creates a shortcut to the transmission_lightcurve dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['transmission_lightcurve']
    except:
        return config_in['transmission']

def from_config_get_transmission(config_in):
    """
    This subroutine creates a shortcut to the transmission_lightcurve dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['transmission']
    except:
        return config_in['transmission_lightcurve']


def from_config_get_system(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    return config_in['system']


def from_config_get_pipeline(config_in):
    """
    This subroutine creates a shortcut to the pipeline parameters dictionary
    :param config_in:
    :return: dictionary
    """
    return config_in['pipeline']


def from_config_get_planet(config_in):
    """
    This subroutine creates a shortcut to the planet dictionary
    :param config_in:
    :return: dictionary
    """
    return config_in['planet']


def from_config_get_star(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    return config_in['star']


def from_config_get_clv_rm(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    return config_in['CLV_RM_correction']


def from_config_get_interstellar_lines(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['interstellar_lines']
    except:
        return None

def from_config_get_molecfit(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['molecfit']
    except:
        return None

def from_config_get_transmission_mcmc(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['transmission_mcmc']
    except:
        return None


def from_config_get_spectral_lines(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['spectral_lines']
    except:
        return None

def from_config_get_interactive_plots(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['interactive_plots']
    except:
        return False

def from_config_get_pca_parameters(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['pca_parameters']
    except:
        return {}


def from_config_get_fullspectrum_parameters(config_in):
    """
    This subroutine creates a shortcut to the system dictionary
    :param config_in:
    :return: dictionary
    """
    try:
        return config_in['full_spectrum']
    except:
        return {}