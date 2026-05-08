from SLOPpy.subroutines.common import *
import os

from SLOPpy.subroutines.mcmc_fit_functions import *
from SLOPpy.subroutines.math_functions import interpolate2d_grid_nocheck

#from SLOPpy.subroutines.interpol import interpolate1d_grid_nocheck

from multiprocessing import Pool
import emcee
import time

# define theta pams
def define_theta_array(model_case,line_iter_dict, planet_dict, n_jitter, allow_emission=False):

    pams_dict = {}  # dictionary containing the index of a given parameter
    pams_list = []  # list wirth the parameter names ordered according to their index
    boundaries = np.empty([0, 2])  # boundaries for MCMC / nested sampling
    theta_start = np.empty(0)  # starting point for MCMC
    lines_center = np.empty(0)  # laboratory wavelength of spectral lines
    pam_index = 0  # keep track of the number of variables

    for line_key, line_val in line_iter_dict['lines'].items():
        pam_name = line_key + '_contrast'
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        if allow_emission:
            boundaries = np.append(boundaries, [[-0.20, 0.20]], axis=0)
        else:
            boundaries = np.append(boundaries, [[0.00, 0.20]], axis=0)
        theta_start = np.append(theta_start, 0.010)
        pam_index += 1

        lines_center = np.append(lines_center, line_val)

        """ skip the inclusion of FWHM as a free parameter for each line
            if the shared FWHM is selected
        """

        if model_case in [0, 1, 2, 3, 10, 11, 14, 20, 21, 24]:
            # if not line_iter_dict['fit_parameters']['shared_fwhm']:
            pam_name = line_key + '_fwhm (km/s)'
            pams_dict[pam_name] = pam_index
            pams_list.append(pam_name)
            boundaries = np.append(boundaries, [[0.00, 100.00]], axis=0)
            theta_start = np.append(theta_start, 5.0)
            pam_index += 1

        # if line_iter_dict['fit_parameters']['fixed_separation']: continue
        # if not line_iter_dict['fit_parameters']['lines_shift']: continue

        if model_case in [0, 2, 10, 12, 20, 22]:
            pam_name = line_key + '_winds  (km/s)'
            pams_dict[pam_name] = pam_index
            pams_list.append(pam_name)
            boundaries = np.append(boundaries, [[-25.00, 25.00]], axis=0)
            theta_start = np.append(theta_start, 0.00)
            pam_index += 1

    if model_case in [12, 13, 15, 22, 23, 25]:
        # if line_iter_dict['fit_parameters']['shared_fwhm']:
        pam_name = 'shared_fwhm  (km/s)'
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        boundaries = np.append(boundaries, [[0.000, 100.00]], axis=0)
        theta_start = np.append(theta_start, 5.000)
        pam_index += 1

    if model_case in [11, 13, 21, 23]:
        # if line_iter_dict['fit_parameters']['fixed_separation'] and line_iter_dict['fit_parameters']['lines_shift']:
        pam_name = 'shared_winds  (km/s)'
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        boundaries = np.append(boundaries, [[-25.0, 25.0]], axis=0)
        theta_start = np.append(theta_start, 0.000)
        pam_index += 1

    if model_case in [0, 1, 10, 11, 12, 13, 14, 15]:
        pams_dict['rp_factor'] = pam_index
        pams_list.append('rp_factor')
        boundaries = np.append(boundaries, [[0.5, 2.0]], axis=0)
        theta_start = np.append(theta_start, 1.0)
        pam_index += 1

    pams_dict['K_planet  (km/s)'] = pam_index
    pams_list.append('K_planet  (km/s)')
    #boundaries = np.append(boundaries,
    #                        [[-300., planet_dict['RV_semiamplitude']
    #                            [0] + 300.]],
    #                        axis=0)
    '''
    boundaries = np.append(boundaries,
                            [[planet_dict['RV_semiamplitude'][0] - 75.,
                              planet_dict['RV_semiamplitude'][0] + 75.]],
                            axis=0)

    '''
    boundaries = np.append(boundaries,
                            [[0.,
                              300.]],
                            axis=0)
    
    theta_start = np.append(
        theta_start, planet_dict['RV_semiamplitude'][0] / 1000.0)
    pam_index += 1


    for i_j in range(0,n_jitter):
        pam_name = 'jitter_' + repr(i_j)
        pams_dict[pam_name] = pam_index
        pams_list.append(pam_name)
        boundaries = np.append(boundaries, [[10**(-12), 0.05]], axis=0)
        theta_start = np.append(theta_start, 10**(-11))
        pam_index += 1

    return lines_center, pams_dict, pams_list, boundaries, theta_start




def emcee_lines_fit_functions(model_case,
                              wave_meshgrid,
                              transmission_spec,
                              transmission_spec_err,
                              clv_rm_radius,
                              clv_rm_grid,
                              planet_RVsinusoid,
                              lines_center,
                              jitter_index,
                              priors_dict,
                              theta_start, boundaries, ndim, nwalkers, ngen, nsteps, nthin, ncpus):

    os.environ["OMP_NUM_THREADS"] = "1"

    #""" Avoid starting values out of boundaries """
    #for nd in range(0, ndim):
    #    sel = (point_star[:,nd] <= boundaries[nd,0]) | (point_star[:,nd] >= boundaries[nd,1])
    #    point_star[sel,nd] = theta_start[nd]


    """
    print(np.shape(wave_append))
    print(np.shape(flux_append))
    print(np.shape(ferr_append))
    print(np.shape(nobs_append))
    print(np.shape(clv_rm_radius))
    print(np.shape(clv_rm_grid))
    print(np.shape(rvsys_PRF2ORF_append))
    print(np.shape(planet_RVsinusoid_append))
    print(np.shape(lines_center))
    """

    model_dictionaries = {
         0: logprob_case00,
         1: logprob_case01,
         2: logprob_case02,
         3: logprob_case03,
        10: logprob_case10,
        11: logprob_case11,
        12: logprob_case12,
        13: logprob_case13,
        14: logprob_case14,
        15: logprob_case15,
        20: logprob_case20,
        21: logprob_case21,
        22: logprob_case22,
        23: logprob_case23,
        24: logprob_case24,
        25: logprob_case25,
    }

    logprob = model_dictionaries[model_case]

    try:
        #from pyde.de import DiffEvol
        from pytransit.utils.de import DiffEvol
        use_pyde = True
    except ImportError:
        print('   Warnign: PyDE is not installed, random initialization point')
        use_pyde = False

    if ngen <= 1 : use_pyde = False

    """ R_p is fixed to 1.0 """
    if model_case in [2, 3, 20, 21, 22, 23, 24, 25]:
        clv_model = interpolate2d_grid_nocheck(1.000, clv_rm_radius, clv_rm_grid)
        args_input = (boundaries,
                        wave_meshgrid,
                        transmission_spec,
                        transmission_spec_err,
                        clv_model,
                        planet_RVsinusoid,
                        lines_center,
                        jitter_index,
                        priors_dict)
    else:
        args_input = (boundaries,
                        wave_meshgrid,
                        transmission_spec,
                        transmission_spec_err,
                        clv_rm_radius,
                        clv_rm_grid,
                        planet_RVsinusoid,
                        lines_center,
                        jitter_index,
                        priors_dict)


    if use_pyde:
        start = time.time()

        with Pool() as pool:

            de = DiffEvol(
                    logprob,
                    boundaries,
                    nwalkers,
                    maximize=True,
                    pool=pool,
                    args=args_input)

            de.optimize(ngen)

            end = time.time()
            print("PyDE global optimization took {0:.1f} seconds".format( end - start))

            theta_start = np.median(de.population, axis=0)
            point_start = de.population
    else:
        point_start = theta_start  + 1e-4 * np.abs(np.random.randn(nwalkers, ndim))
        point_start[0, :] = theta_start

    start = time.time()
    with Pool(processes=ncpus) as pool:

        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        logprob,
                                        args=args_input,
                                        pool=pool)

        population, prob, state = sampler.run_mcmc(point_start,
                                                   nsteps,
                                                   thin=nthin,
                                                   progress=True)

    end = time.time()
    print()
    print("emcee MCMC optimization took {0:.1f} seconds".format(end - start))

    return population, sampler.chain, sampler.lnprobability, point_start


def return_model(model_case,
                theta,
                wave_meshgrid,
                clv_rm_radius,
                clv_rm_grid,
                planet_RVsinusoid,
                lines_center,
                jitter_index):

    ndim = len(theta)
    boundaries = np.empty([ndim, 2])
    boundaries[:,0] = theta - 1.
    boundaries[:,1] = theta + 1.
    transmission_spec = np.ones(np.shape(wave_meshgrid))
    transmission_spec_err = np.ones(np.shape(wave_meshgrid))

    model_dictionaries = {
         0: logprob_case00,
         1: logprob_case01,
         2: logprob_case02,
         3: logprob_case03,
        10: logprob_case10,
        11: logprob_case11,
        12: logprob_case12,
        13: logprob_case13,
        14: logprob_case14,
        15: logprob_case15,
        20: logprob_case20,
        21: logprob_case21,
        22: logprob_case22,
        23: logprob_case23,
        24: logprob_case24,
        25: logprob_case25,
    }

    logprob = model_dictionaries[model_case]

    if model_case in [2, 3, 20, 21, 22, 23, 24, 25]:
        clv_model = interpolate2d_grid_nocheck(1.000, clv_rm_radius, clv_rm_grid)
        lines_model, _, lines_array, planet_K, planet_R, jitter = logprob(
                theta,
                boundaries,
                wave_meshgrid,
                transmission_spec,
                transmission_spec_err,
                clv_model,
                planet_RVsinusoid,
                lines_center,
                jitter_index,
                {},
                return_models=True)

    else:
        lines_model, clv_model, lines_array, planet_K, planet_R, jitter = logprob(
                theta,
                boundaries,
                wave_meshgrid,
                transmission_spec,
                transmission_spec_err,
                clv_rm_radius,
                clv_rm_grid,
                planet_RVsinusoid,
                lines_center,
                jitter_index,
                {},
                return_models=True)


    return lines_model, clv_model, lines_array, planet_K, planet_R, jitter


def emcee_flatten_median(population, sampler_chain, sampler_lnprobability,  nburnin, nthin, nwalkers):

    flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
    flat_lnprob, _ = emcee_flatlnprob(
            sampler_lnprobability, nburnin, nthin, population, nwalkers)

    lnprob_med = compute_value_sigma(flat_lnprob)
    chain_med = compute_value_sigma(flat_chain)
    chain_MAP, lnprob_MAP = pick_MAP_parameters(flat_chain, flat_lnprob)

    #n_samplings, n_pams = np.shape(flat_chain)
    return flat_chain, flat_lnprob, chain_med, chain_MAP, lnprob_med, lnprob_MAP

def emcee_compute_BIC_AIC(lnprob_med, lnprob_MAP, ndata, ndim):

    print()
    print(' LN posterior: {0:12f}   {1:12f} {2:12f} (15-84 p)  MAP:  {0:12f}'.format(
        lnprob_med[0], lnprob_med[2], lnprob_med[1], lnprob_MAP))

    BIC = -2.0 * lnprob_med[0] + np.log(ndata) * ndim
    AIC = -2.0 * lnprob_med[0] + 2.0 * ndim
    AICc = AIC + (2.0 + 2.0 * ndim) * ndim / (ndata - ndim - 1.0)

    print()
    print(' Median BIC    = {}'.format(BIC))
    print(' Median AIC    = {}'.format(AIC))
    print(' Median AICc   = {}'.format(AICc))

    BIC_map = -2.0 * lnprob_MAP + np.log(ndata) * ndim
    AIC_map = -2.0 * lnprob_MAP + 2.0 * ndim
    AICc_map = AIC + (2.0 + 2.0 * ndim) * ndim / (ndata - ndim - 1.0)

    print()
    print(' MAP BIC       = {}'.format(BIC_map))
    print(' MAP AIC       = {}'.format(AIC_map))
    print(' MAP AICc      = {}'.format(AICc_map))


def emcee_burnin_check(chain, nburnin, nthin, nwalkers=False):
    nburn = int(nburnin / nthin)
    modified = False

    if not nwalkers:
        _, d, _ = np.shape(chain)
    else:
        v1, v2 = np.shape(chain)
        if v1 == nwalkers:
            d = v2
        else:
            d = v1

    if nburn >= d * 0.9:
        nburn = int(d / 4)
        modified = True

    return nburn, modified

def emcee_flatchain(chain, nburnin, nthin):
    """flattening of the emcee chains with removal of burn-in"""
    nburn, _ = emcee_burnin_check(chain, nburnin, nthin)
    s = chain[:, nburn:, :].shape
    return chain[:, nburn:, :].reshape(s[0] * s[1], s[2])


def emcee_flatlnprob(lnprob, nburnin, nthin, population, nwalkers):

    nburn, _  = emcee_burnin_check(lnprob, nburnin, nthin, nwalkers)

    v1, v2 = np.shape(lnprob)
    if v1 == nwalkers:
        s = lnprob[:, nburn:].shape
        return lnprob[:, nburn:].reshape(s[0] * s[1]), lnprob.T
    else:
        s = lnprob[nburn:, :].shape
        return lnprob[nburn:, :].reshape(s[0] * s[1]), lnprob

def GelmanRubin_v2(sampler_chain):
    """
    :param chain_T:
    :return:
    """

    """
    from http://joergdietrich.github.io/emcee-convergence.html
    """
    ssq = np.var(sampler_chain, axis=1, ddof=1)
    W = np.mean(ssq, axis=0)
    theta_b = np.mean(sampler_chain, axis=1)
    theta_bb = np.mean(theta_b, axis=0)
    m = sampler_chain.shape[0] * 1.0
    n = sampler_chain.shape[1] * 1.0
    B = n / (m - 1) * np.sum((theta_bb - theta_b) ** 2, axis=0)
    var_theta = (n - 1) / n * W + 1 / n * B
    Rhat = np.sqrt(var_theta / W)
    return Rhat


def compute_value_sigma(samples):
    if np.size(np.shape(samples)) == 1:
        sample_med = np.zeros(3)
        finite_samples = samples[np.isfinite(samples)]
        if len(finite_samples) == 0:
            sample_med[:] = np.nan
            return sample_med
        sample_tmp = np.percentile(finite_samples, [15.865, 50, 84.135], axis=0)

        sample_med[0] = sample_tmp[1]
        sample_med[1] = sample_tmp[2] - sample_tmp[1]
        sample_med[2] = sample_tmp[1] - sample_tmp[0]

    elif np.size(np.shape(samples)) == 2:
        #sample_med = np.asarray(list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                   #zip(*np.percentile(samples, [15.865, 50, 84.135], axis=0)))))
        sample_med = np.zeros((samples.shape[1],3))
        for k in range(samples.shape[1]):
            ttt = samples[:,k]
            sample_tmp = np.percentile(ttt[np.isfinite(ttt)], [15.865, 50, 84.135], axis=0)
            sample_med[k,0] = sample_tmp[1]
            sample_med[k,1] = sample_tmp[2] - sample_tmp[1]
            sample_med[k,2] = sample_tmp[1] - sample_tmp[0]


    else:
        print('ERROR!!! ')
        return None
    return sample_med

def pick_MAP_parameters(samples, lnprob):

    indmax = np.argmax(lnprob)
    if np.size(np.shape(samples)) == 1:
        return samples[indmax], lnprob[indmax]
    elif np.size(np.shape(samples)) == 2:
        return samples[indmax, :], lnprob[indmax]
    else:
        print('ERROR!!! ')
        return None
