from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.constants import *
from SLOPpy.subroutines.kepler_exo import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.math_functions import *

from astropy.convolution import Gaussian1DKernel, convolve

__all__ = ['compute_clv_rm_models', 'plot_clv_rm_models']

subroutine_name = 'clv_rm_models'


def compute_clv_rm_models(config_in):
    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)

    planet_dict = from_config_get_planet(config_in)
    star_dict = from_config_get_star(config_in)
    clv_rm_dict = from_config_get_clv_rm(config_in)

    spectral_lines = from_config_get_spectral_lines(config_in)

    # un-convolved portion of the spectrum given by range_boundaries +-
    wave_fix_convo = 1.0

    # Added back-compatibility to old or "wrong" keys
    norm_dict = clv_rm_dict.get('normalization', {})
    norm_pams={}
    norm_pams['model_poly_degree'] = norm_dict.get('model_poly_degree', 2)
    norm_pams['spectra_poly_degree'] = norm_dict.get('spectra_poly_degree', 2)
    norm_pams['lower_threshold'] = norm_dict.get('lower_threshold', 0.950)
    norm_pams['percentile_selection'] = norm_dict.get('percentile_selection', 10)

    try:
        synthesis = load_from_cpickle('clv_rm_synthesis', config_in['output'])
        star_grid = load_from_cpickle('clv_rm_star_grid', config_in['output'])

        if not config_in['settings'].get('full_output', False):
            for night in night_dict:
                clv_rm_models = load_from_cpickle(
                    'clv_rm_models', config_in['output'], night)

        print("{0:45s}                         {1:s}".format(
            subroutine_name, 'Retrieved'))


    except:
        print("{0:45s}                         {1:s}".format(
            subroutine_name, 'Computing'))
        print()

        """
        Loading the spectral synthesis results, at the moment only SME output is supported.
        Properties of the synthesis data files
        - limb_angles: this is an input to SME, so it is specific on how the synthesis has been performed
        - spectra: stellar spectrum as a function of the limb angle, sampled near the spectral lines
        - model: integrated spectrum of the star
        """

        synthesis_data_limb_angles = np.genfromtxt(
            clv_rm_dict['synthesis_files'] + '_muvals.txt', dtype=np.single)
        synthesis_data_spectra = np.genfromtxt(
            clv_rm_dict['synthesis_files'] + '_spectra.txt', dtype=np.double)
        synthesis_data_model = np.genfromtxt(
            clv_rm_dict['synthesis_files'] + '_model.txt', dtype=np.double)

        synthesis = {
            'surface': {
                'wave': synthesis_data_spectra[:, 0],
                'flux': np.single(synthesis_data_spectra[:, 1:]),
                'n_mu': np.size(synthesis_data_limb_angles),
                'mu': synthesis_data_limb_angles
            },
            'total': {
                'wave': synthesis_data_model[:, 0],
                'norm': np.single(synthesis_data_model[:, 1]),
            }
        }

        """ Setting up the array for model computation """

        synthesis['total']['step'] = synthesis['total']['wave'] * 0.0
        synthesis['total']['step'][1:] = synthesis['total']['wave'][1:] - \
            synthesis['total']['wave'][:-1]
        synthesis['total']['step'][0] = synthesis['total']['step'][1]

        synthesis['surface']['step'] = synthesis['surface']['wave'] * 0.0
        synthesis['surface']['step'][1:] = synthesis['surface']['wave'][1:] - \
            synthesis['surface']['wave'][:-1]
        synthesis['surface']['step'][0] = synthesis['surface']['step'][1]

        synthesis['surface']['wave_out'] = np.arange(synthesis['surface']['wave'][0],
                                                     synthesis['surface']['wave'][-1],
                                                     clv_rm_dict['rebinning_step'])
        synthesis['surface']['size_out'] = np.size(
            synthesis['surface']['wave_out'], axis=0)
        synthesis['surface']['step_out'] = np.ones(
            synthesis['surface']['size_out']) * clv_rm_dict['rebinning_step']

        synthesis['total']['norm_out'] = rebin_1d_to_1d(synthesis['total']['wave'],
                                                        synthesis['total']['step'],
                                                        synthesis['total']['norm'],
                                                        synthesis['surface']['wave_out'],
                                                        synthesis['surface']['step_out'],
                                                        method='exact_flux',
                                                        preserve_flux=False)

        """  Check if the number of spectra corresponds to the number of limb angle values """
        if np.size(synthesis['surface']['flux'], axis=1) != synthesis['surface']['n_mu']:
            print('ERROR in loading the stellar spectra')

        """
        Setting up the grid of stellar spectra for the CLV and RM computation
        odd number of points to include the zero value
        """
        star_grid = {
            'n_grid': clv_rm_dict['n_gridpoints'],
            'half_grid': int((clv_rm_dict['n_gridpoints'] - 1) / 2)
        }

        """ Coordinates of the centers of each grid cell (add offset) """
        star_grid['xx'] = np.linspace(-1.0000000000000, 1.0000000000000,
                                      star_grid['n_grid'], dtype=np.double)
        star_grid['xc'], star_grid['yc'] = np.meshgrid(
            star_grid['xx'], star_grid['xx'], indexing='xy')
        # check the Note section of the wiki page of meshgrid
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html

        """ Distance of each grid cell from the center of the stellar disk """
        star_grid['rc'] = np.sqrt(star_grid['xc'] ** 2 + star_grid['yc'] ** 2)
        # Must avoid negative numbers inside the square root
        star_grid['inside'] = star_grid['rc'] < 1.0000000000000
        # Must avoid negative numbers inside the square root
        star_grid['outside'] = star_grid['rc'] >= 1.00000000000000

        """ Determine the mu angle for each grid cell, as a function of radius. """
        star_grid['mu'] = np.zeros([star_grid['n_grid'], star_grid['n_grid']],
                                   dtype=np.double)  # initialization of the matrix with the mu values
        star_grid['mu'][star_grid['inside']] = np.sqrt(
            1. - star_grid['rc'][star_grid['inside']] ** 2)

        """  2.2 Determine the Doppler shift to apply to the spectrum of each grid cell, from Cegla+2016 """

        star_grid['x_ortho'] = star_grid['xc'] * np.cos(star_dict['lambda'][0] * deg2rad) \
            - star_grid['yc'] * np.sin(
            star_dict['lambda'][0] * deg2rad)  # orthogonal distances from the spin-axis
        star_grid['y_ortho'] = star_grid['xc'] * np.sin(star_dict['lambda'][0] * deg2rad) \
            + star_grid['yc'] * np.cos(star_dict['lambda'][0] * deg2rad)

        star_grid['r_ortho'] = np.sqrt(
            star_grid['x_ortho'] ** 2 + star_grid['y_ortho'] ** 2)
        star_grid['z_ortho'] = np.zeros([star_grid['n_grid'], star_grid['n_grid']],
                                        dtype=np.double)  # initialization of the matrix
        star_grid['z_ortho'][star_grid['inside']] = np.sqrt(
            1. -star_grid['r_ortho'][star_grid['inside']] ** 2)

        """ rotate the coordinate system around the x_ortho axis by an angle: """
        star_grid['beta'] = (np.pi / 2.) - \
            star_dict['inclination'][0] * deg2rad

        """ orthogonal distance from the stellar equator """
        star_grid['yp_ortho'] = star_grid['z_ortho'] * np.sin(star_grid['beta']) + star_grid['y_ortho'] * np.cos(
            star_grid['beta'])

        """ stellar rotational velocity for a given position """
        star_grid['v_star'] = star_grid['x_ortho'] * star_dict['vsini'][0] * (
            1. -star_dict['alpha'][0] * star_grid['yp_ortho'] ** 2)
        # Null velocity for points outside the stellar surface
        star_grid['v_star'][star_grid['outside']] = 0.0

        """ Associate a synthetic spectrum to each cell """

    """ recomputation of spectra_mu  - most likely it has been deleted from the
    output file
    """

    star_grid['spectra_mu'] = [[0] * star_grid['n_grid']
                                for i in range(star_grid['n_grid'])]

    for x in range(0, star_grid['n_grid']):
        for y in range(0, star_grid['n_grid']):

            if star_grid['outside'][y, x]:
                continue

            index_closer = np.abs(
                synthesis['surface']['mu'] - star_grid['mu'][y, x]).argmin()  # take the index of the closer value

            if star_grid['mu'][y, x] in synthesis['surface']['mu']:
                star_grid['spectra_mu'][x][y] = synthesis['surface']['flux'][:, index_closer]
                continue
            elif index_closer == synthesis['surface']['n_mu'] - 1 or \
                    synthesis['surface']['mu'][index_closer] > star_grid['mu'][y, x]:
                mu_ind0 = index_closer - 1
                mu_ind1 = index_closer
            else:
                mu_ind0 = index_closer
                mu_ind1 = index_closer + 1

            diff_mu = synthesis['surface']['mu'][mu_ind1] - \
                synthesis['surface']['mu'][mu_ind0]

            star_grid['spectra_mu'][x][y] = synthesis['surface']['flux'][:, mu_ind0] \
                + (star_grid['mu'][y, x] - synthesis['surface']['mu'][mu_ind0]) / diff_mu \
                * (synthesis['surface']['flux'][:, mu_ind1]
                   - synthesis['surface']['flux'][:, mu_ind0])

    """ Computation of the continuum level (total flux is already normalized)"""
    star_grid['continuum'] = [[0] * star_grid['n_grid']
                              for i in range(star_grid['n_grid'])]

    spectral_window = ((synthesis['surface']['wave'] > clv_rm_dict['continuum_range'][0]) &
                       (synthesis['surface']['wave'] < clv_rm_dict['continuum_range'][1]))

    for x in range(0, star_grid['n_grid']):
        for y in range(0, star_grid['n_grid']):

            if star_grid['outside'][y, x]:
                continue
            star_grid['continuum'][x][y] = np.median(
                star_grid['spectra_mu'][x][y][spectral_window])

    star_grid['continuum_level'] = np.sum(star_grid['continuum'])

    """
    Setting up the grid for the rescaling factor of the planetary radius
    """
    try:
        radius_grid = np.arange(clv_rm_dict['radius_factor'][0],
                                clv_rm_dict['radius_factor'][1] +
                                clv_rm_dict['radius_factor'][2],
                                clv_rm_dict['radius_factor'][2])
    except KeyError:
        radius_grid = np.arange(0.5, 2.6, 0.1)

    for night in night_dict:
        """ Retrieving the list of observations"""

        print()
        print('compute_CLV_RM_models                   Night: ', night)
        try:
            clv_rm_models = load_from_cpickle(
                'clv_rm_models', config_in['output'], night)
            continue
        except:
            print()
            print('  No CLV & RM correction files found, computing now ')

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle(
            'observational_pams', config_in['output'], night)

        instrument = night_dict[night]['instrument']

        clv_rm_models = {
            'common': {
                'wave': synthesis['surface']['wave_out'],
                'step': synthesis['surface']['step_out'],
                'norm': synthesis['total']['norm_out'],
                'continuum_level': star_grid['continuum_level'],
                'radius_grid': radius_grid,
                'n_radius_grid': len(radius_grid)
            }
        }

        clv_rm_models['common']['convolution_dlambda'] = \
            np.median(clv_rm_models['common']['wave']) / \
            instrument_dict[instrument]['resolution']
        clv_rm_models['common']['convolution_sigma'] = \
            clv_rm_models['common']['convolution_dlambda'] / \
            np.median(clv_rm_models['common']['step'])
        gaussian = Gaussian1DKernel(
            stddev=clv_rm_models['common']['convolution_sigma'])

        clv_rm_models['common']['norm_convolved'] = convolve(
            clv_rm_models['common']['norm'], gaussian)

        """ Fixing border effect (we took already wave_extension angstrom outside of the
        actual range, so doing it this way is fine)"""

        wave_fix_convolution = (clv_rm_models['common']['wave'] > clv_rm_models['common']['wave'][0]+wave_fix_convo) \
            | (clv_rm_models['common']['wave'] > clv_rm_models['common']['wave'][-1]-wave_fix_convo)
        clv_rm_models['common']['norm_convolved'][wave_fix_convolution] = clv_rm_models['common']['norm'][wave_fix_convolution]

        """
        Computation of the first derivative, useful to identify
        continuum level. This method is prone to errors for
        observational data, but it's quite robust for synthetic spectra
        if jumps in wavelngth are small
        """

        clv_rm_models['common']['norm_convolved_derivative'] = \
                first_derivative(clv_rm_models['common']['wave'],
                                    clv_rm_models['common']['norm_convolved'])


        # Using only the 10percentile of values of the derivative around zero
        cont_10perc = np.percentile(np.abs(clv_rm_models['common']['norm_convolved_derivative']), norm_pams['percentile_selection'])

        clv_rm_models['common']['norm_convolved_bool'] = (np.abs(clv_rm_models['common']['norm_convolved_derivative']) < cont_10perc) \
            & (clv_rm_models['common']['norm_convolved']> norm_pams['lower_threshold'])

        print('  Number of points within 10percentile: {0:10.0f}'.format(np.sum((np.abs(clv_rm_models['common']['norm_convolved_derivative']) < cont_10perc))))
        print('  Number of points above threshold: {0:10.0f}'.format(np.sum( (clv_rm_models['common']['norm_convolved']> norm_pams['lower_threshold']))))

        norm_convolved_bool = (np.abs(clv_rm_models['common']['norm_convolved_derivative']) < cont_10perc) \
            & (clv_rm_models['common']['norm_convolved']> norm_pams['lower_threshold'])

        if np.sum(norm_convolved_bool) < 100:
            print('  Lower threshold decreased by 80% to allow point selection ', norm_pams['lower_threshold']*0.80)

            clv_rm_models['common']['norm_convolved_bool'] = (np.abs(clv_rm_models['common']['norm_convolved_derivative']) < cont_10perc) \
                & (clv_rm_models['common']['norm_convolved']> norm_pams['lower_threshold']*0.80)
        else:
            clv_rm_models['common']['norm_convolved_bool'] = norm_convolved_bool


        processed = {}

        print()
        for obs in lists['observations']:

            print('  Computing CLV+RM correction for ', obs)

            processed[obs] = {}
            clv_rm_models[obs] = {}

            n_oversampling = int(
                observational_pams[obs]['EXPTIME'] / clv_rm_dict['time_step'])

            if n_oversampling % 2 == 0:
                n_oversampling += 1
            half_time = observational_pams[obs]['EXPTIME'] / 2 / 86400.

            processed[obs]['bjd_oversampling'] = np.linspace(observational_pams[obs]['BJD'] - half_time,
                                                             observational_pams[obs]['BJD'] + half_time,
                                                             n_oversampling, dtype=np.double)
            if planet_dict['orbit'] == 'circular':
                # Time of pericenter concides with transit time, if we assume e=0 and omega=np.pi/2.
                eccentricity = 0.00
                omega_rad = np.pi / 2.
                # Tcent is assumed as reference time
                Tref = planet_dict['reference_time_of_transit'][0]
                Tcent_Tref = 0.000
            else:
                omega_rad = planet_dict['omega'][0] * deg2rad
                Tref = planet_dict.get('reference_time', planet_dict['reference_time_of_transit'][0])
                Tcent_Tref = planet_dict['reference_time_of_transit'][0] - Tref
                eccentricity = planet_dict['eccentricity'][0]
            inclination_rad = planet_dict['inclination'][0] * deg2rad

            true_anomaly, orbital_distance_ratio = kepler_true_anomaly_orbital_distance(
                processed[obs]['bjd_oversampling'] - Tref,
                Tcent_Tref,
                planet_dict['period'][0],
                eccentricity,
                omega_rad,
                planet_dict['semimajor_axis_ratio'][0])

            """ planet position during its orbital motion, in unit of stellar radius"""
            # Following Murray+Correia 2011 , with the argument of the ascending node set to zero.
            # 1) the ascending node coincide with the X axis
            # 2) the reference plance coincide with the plane of the sky
            processed[obs]['planet_position'] = {
                'xp': -orbital_distance_ratio * (np.cos(omega_rad + true_anomaly)),
                'yp': orbital_distance_ratio * (np.sin(omega_rad + true_anomaly) * np.cos(inclination_rad)),
                'zp': orbital_distance_ratio * (np.sin(inclination_rad) * np.sin(omega_rad + true_anomaly))
            }
            # projected distance of the planet's center to the stellar center
            processed[obs]['planet_position']['rp'] = np.sqrt(processed[obs]['planet_position']['xp'] ** 2
                                                              + processed[obs]['planet_position']['yp'] ** 2)

            # obscured flux integrated over the full epoch
            # grid n_radius_grid X size_out (of spectral model)
            clv_rm_models[obs]['missing_flux'] = np.zeros(
                [len(radius_grid), synthesis['surface']['size_out']], dtype=np.double)

            # iterating on the sub-exposures
            for j, zeta in enumerate(processed[obs]['planet_position']['zp']):

                if zeta > 0 and processed[obs]['planet_position']['rp'][j] < 1. + planet_dict['radius_ratio'][0]:
                    # the planet is in the foreground or inside the stellar disk, continue
                    # adjustment: computation is performed even if only part of the planet is shadowing the star

                    rd = np.sqrt((processed[obs]['planet_position']['xp'][j] - star_grid['xc']) ** 2 +
                                 (processed[obs]['planet_position']['yp'][j] - star_grid['yc']) ** 2)

                    # iterating on the cell grid
                    for x in range(0, star_grid['n_grid']):
                        for y in range(0, star_grid['n_grid']):

                            # skip the step if the cell is outside the stellar disk
                            # or if the cell is not shadowed by the planet when the largest possible size is considered
                            if star_grid['outside'][y, x] or rd[y, x] > planet_dict['radius_ratio'][0]*radius_grid[-1]:
                                continue

                            # rescaled planetary radius selection
                            grid_sel = (
                                rd[y, x] <= planet_dict['radius_ratio'][0]*radius_grid)
                            # stellar flux in the masked region
                            flux_tmp = rebin_1d_to_1d(synthesis['surface']['wave'],
                                                      synthesis['surface']['step'],
                                                      star_grid['spectra_mu'][x][y],
                                                      clv_rm_models['common']['wave'],
                                                      clv_rm_models['common']['step'],
                                                      rv_shift=star_grid['v_star'][y, x],
                                                      method='exact_flux',
                                                      preserve_flux=False)

                            # fixing zero values that may have been introduced by
                            # the rebinning process from an extremely irregular sampling

                            ind_sel = np.where(flux_tmp < 0.)[0]
                            for ii in ind_sel:
                                if ii == 0:
                                    flux_tmp[ii] = flux_tmp[ii + 1]
                                elif ii == np.size(flux_tmp) - 1:
                                    flux_tmp[ii] = flux_tmp[ii - 1]
                                else:
                                    flux_tmp[ii] = (
                                        flux_tmp[ii - 1] + flux_tmp[ii + 1]) / 2.

                            """
                            Outer product of the radius selection array (size=M)
                            and the flux array (N) so that it can be summed
                            properly to the MxN missing_flux matrix.
                            """
                            clv_rm_models[obs]['missing_flux'] += \
                                np.outer(grid_sel, flux_tmp)

            clv_rm_models[obs]['missing_flux'] /= n_oversampling
            clv_rm_models[obs]['stellar_spectra'] = \
                np.outer(np.ones(len(radius_grid)), clv_rm_models['common']['norm']) \
                - (clv_rm_models[obs]['missing_flux'] /
                   clv_rm_models['common']['continuum_level'])

            clv_rm_models[obs]['stellar_spectra_convolved'] = \
                np.zeros([len(radius_grid), synthesis['surface']['size_out']],
                         dtype=np.double)

            clv_rm_models[obs]['clv_rm_model_convolved'] = \
                np.zeros([len(radius_grid), synthesis['surface']['size_out']],
                         dtype=np.double)

            clv_rm_models[obs]['clv_rm_model_convolved_derivative'] = \
                np.zeros([len(radius_grid), synthesis['surface']['size_out']],
                         dtype=np.double)

            clv_rm_models[obs]['clv_rm_model_convolved_continuum_bool'] = \
                np.zeros([len(radius_grid), synthesis['surface']['size_out']],
                         dtype=bool)

            clv_rm_models[obs]['clv_rm_model_convolved_normalized'] = \
                np.zeros([len(radius_grid), synthesis['surface']['size_out']],
                         dtype=np.double)



            for ii in range(0, len(radius_grid)):
                clv_rm_models[obs]['stellar_spectra_convolved'][ii, :] = \
                    convolve(clv_rm_models[obs]['stellar_spectra'][ii, :],
                             gaussian)

                clv_rm_models[obs]['stellar_spectra_convolved'][ii, wave_fix_convolution] = \
                    clv_rm_models[obs]['stellar_spectra'][ii, wave_fix_convolution]

                """
                This is the theoretical transmission spectrum in the stellar reference frame
                when only CLV and RM effects are present (no atmospheric
                transmission)
                """
                clv_rm_models[obs]['clv_rm_model_convolved'][ii, :] = \
                    clv_rm_models[obs]['stellar_spectra_convolved'][ii, :] \
                    / clv_rm_models['common']['norm_convolved']

                """
                High-resolution transmission spectra are always rescaled for
                their continuum because in fiber-fed spectrographs the
                information on the absolute flux of the star is lost.

                If not using the normalized spectrum, normalization factor must
                be included somehow when correcting for the CLV+RM, before
                fitting the atomic absoprtion lines
                """

                normalization_function = np.polynomial.chebyshev.Chebyshev.fit(
                    clv_rm_models['common']['wave'][clv_rm_models['common']['norm_convolved_bool']],
                    clv_rm_models[obs]['clv_rm_model_convolved'][ii, :][clv_rm_models['common']['norm_convolved_bool']],
                    deg=norm_pams['model_poly_degree']
                )

                clv_rm_models[obs]['clv_rm_model_convolved_normalized'][ii, :] = clv_rm_models[obs]['clv_rm_model_convolved'][ii, :] / normalization_function(clv_rm_models['common']['wave'])

                #plt.plot(clv_rm_models['common']['wave'], clv_rm_models[obs]['clv_rm_model_convolved_normalized'][ii, :])
                #plt.plot(clv_rm_models['common']['wave'], clv_rm_models[obs]['clv_rm_model_convolved'][ii, :])

            #plt.show()

            """ In the planetary reference frame, the corrected transmission
            spectrum T_corr is given by

            T_corr = T_input * (synthetic_convolved /
            stellar_spectra_convolved),

            where: T_input: transmission spectrum before the correction
            synthetic_convolved: integrated synthetic stellar spectrum,
            convolved for the instrumental resolution.
            stellar_spectra_convolved: stellar spectrum after removing the
            contribute of the stellar surface covered by the planet, convolved
            for the instrumental resolution (synthetic_convolved and
            stellar_spectra_convolved are in the stellar rest frame must be
            rebinned in the planetary rest frame)

            Since clv_rm_model_convolved = stellar_spectra_convolved /
            synthetic_convolved the observed transmission spectrum must be
            DIVIDED by clv_rm_model_convolved
            """


        save_to_cpickle('clv_rm_models', clv_rm_models,
                        config_in['output'], night)
        clv_rm_models = None  # Forcing memory de-allocation


    if not config_in['settings'].get('full_output', False):
        del star_grid['spectra_mu']

    save_to_cpickle('clv_rm_star_grid', star_grid, config_in['output'])
    save_to_cpickle('clv_rm_synthesis', synthesis, config_in['output'])


def plot_clv_rm_models(config_in, night_input=''):
    night_dict = from_config_get_nights(config_in)

    synthesis = load_from_cpickle('clv_rm_synthesis', config_in['output'])
    star_grid = load_from_cpickle('clv_rm_star_grid', config_in['output'])

    if night_input == '':
        # Visualize the mu of star
        fig = plt.figure(figsize=(8, 6.5))

        plt.title('Limb angle')
        plt.contourf(star_grid['xx'], star_grid['xx'],
                     star_grid['mu'], 60, cmap=plt.cm.viridis)
        plt.colorbar(label='$\mu$')  # draw colorbar
        # plot data points.
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.xlabel('x [R_s]')
        plt.ylabel('y [R_s]')
        plt.show()

        # Visualize the RV of star
        fig = plt.figure(figsize=(8, 6.5))

        # CS = plt.contour(xx,xx,v_star,50,linewidths=0.5,colors='k')
        plt.title('Radial velocity field')

        plt.contourf(star_grid['xx'], star_grid['xx'],
                     star_grid['v_star'], 100, cmap=plt.cm.seismic)
        plt.colorbar(label='v_star')  # draw colorbar
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.xlabel('x [R_s]')
        plt.ylabel('y [R_s]')
        plt.show()

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        lists = load_from_cpickle('lists', config_in['output'], night)
        clv_rm_models = load_from_cpickle(
            'clv_rm_models', config_in['output'], night)
        observational_pams = load_from_cpickle(
            'observational_pams', config_in['output'], night)

        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(
            lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        i0_radius = np.argmin(
            np.abs(clv_rm_models['common']['radius_grid']-1.00))

        for obs in lists['transit_in']:
            ax1.plot(clv_rm_models['common']['wave'],
                     clv_rm_models[obs]['stellar_spectra'][i0_radius, :],
                     color=colors_plot['mBJD'][obs], alpha=0.2)

            ax1.plot(clv_rm_models['common']['wave'],
                     clv_rm_models[obs]['missing_flux'][i0_radius, :] /
                     clv_rm_models['common']['continuum_level'],
                     color=colors_plot['mBJD'][obs], alpha=0.2)

            ax2.plot(clv_rm_models['common']['wave'],
                     clv_rm_models[obs]['stellar_spectra'][-1, :],
                     color=colors_plot['mBJD'][obs], alpha=0.2)

            ax2.plot(clv_rm_models['common']['wave'],
                     clv_rm_models[obs]['missing_flux'][-1, :] /
                     clv_rm_models['common']['continuum_level'],
                     color=colors_plot['mBJD'][obs], alpha=0.2)

        # for obs in lists['transit_out']:
        #    ax2.plot(clv_rm_models['common']['wave'],
        #             clv_rm_models[obs]['stellar_spectra'],
        #             color=colors_plot['mBJD'][obs], alpha=0.2)

        ax1.set_title(
            'Night: {0:s} \n Input spectra, stellar radius'.format(night))
        ax2.set_title('Stellar radius x {0:2.2f}'.format(
            clv_rm_models['common']['radius_grid'][-1]))
        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(
            cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
