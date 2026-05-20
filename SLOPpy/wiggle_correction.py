from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.kepler_exo import compute_planet_RV

from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline
import glob
import os

__all__ = ['compute_wiggle_correction',
           'plot_wiggle_correction']

subroutine_name = 'wiggle_correction'


def _build_mask(wave, mask_regions):
    """Return a boolean array: True where wavelength is NOT masked (good for fitting)."""
    good = np.ones(len(wave), dtype=bool)
    for region in mask_regions:
        good &= ~((wave >= region[0]) & (wave <= region[1]))
    return good


def _spline_wiggle_correction(wave, flux, flux_err, mask_regions,
                              smoothing_factor=None, degree=3,
                              knot_spacing=2.0):
    """
    Correct wiggles in a 1-D spectrum by fitting a spline to the continuum
    (spectral-line regions masked) and dividing by it.

    Two fitting modes, selected by ``knot_spacing``:

    **Fixed-knot mode** (recommended, ``knot_spacing > 0``):
        Uses ``LSQUnivariateSpline`` with interior knots placed every
        ``knot_spacing`` Angstrom.  No weights or smoothing parameter needed.
        The spline captures features whose period > 2 × knot_spacing Å while
        ignoring narrow lines and pixel-to-pixel noise.  This is robust
        regardless of the noise estimate stored in ``flux_err``.

    **Auto-smoothing mode** (``knot_spacing=None`` or 0):
        Falls back to ``UnivariateSpline`` with ``s=smoothing_factor``.
        *Not recommended* for ESPRESSO: the propagated ``deblazed_err`` is
        comparable to the wiggle amplitude, causing scipy to fit a constant.

    Parameters
    ----------
    wave : array_like
        Wavelength array (Angstrom).
    flux : array_like
        Flux (or transmission ratio) array, normalised around 1.
    flux_err : array_like or None
        Flux error array (used only in auto-smoothing mode).
    mask_regions : list of [float, float]
        Wavelength intervals [wave_min, wave_max] to exclude from the fit.
    smoothing_factor : float or None
        Passed to ``UnivariateSpline`` only in auto-smoothing mode.
    degree : int
        Spline degree (default 3, cubic).
    knot_spacing : float or None
        Spacing (Å) between interior knots.  Default 2.0 Å.
        Set to ``None`` or 0 to fall back to auto-smoothing mode.

    Returns
    -------
    corrected : ndarray  – flux divided by the spline continuum
    corrected_err : ndarray or None
    continuum : ndarray  – the fitted spline evaluated on the full wavelength
    """
    wave = np.asarray(wave, dtype=np.float64)
    flux = np.asarray(flux, dtype=np.float64)

    good = _build_mask(wave, mask_regions)
    finite = np.isfinite(flux) & (flux > 0)
    sel = good & finite

    if np.sum(sel) < degree + 1:
        continuum = np.ones_like(flux)
        if flux_err is not None:
            return flux.copy(), np.asarray(flux_err, dtype=np.float64).copy(), continuum
        return flux.copy(), None, continuum

    use_fixed_knots = knot_spacing is not None and float(knot_spacing) > 0

    if use_fixed_knots:
        # ---------------------------------------------------------------
        # Fixed-knot LSQ spline – robust against noise-scale ambiguity
        # ---------------------------------------------------------------
        knot_spacing = float(knot_spacing)
        w_min = wave[sel].min()
        w_max = wave[sel].max()
        span = w_max - w_min

        # Number of interior knots: span / knot_spacing, but at least degree
        n_interior = max(degree, int(np.floor(span / knot_spacing)) - 1)
        # Place knots strictly inside the data range
        t_interior = np.linspace(w_min + knot_spacing * 0.5,
                                 w_max - knot_spacing * 0.5,
                                 n_interior)
        # Keep only knots strictly inside data bounds
        t_interior = t_interior[(t_interior > w_min) & (t_interior < w_max)]

        if len(t_interior) < 1:
            # Span too small for even one interior knot – fall back to polynomial
            t_interior = np.array([(w_min + w_max) / 2.0])

        try:
            spl = LSQUnivariateSpline(wave[sel], flux[sel],
                                      t=t_interior, k=degree)
        except Exception:
            # Any LSQ failure → return unchanged
            continuum = np.ones_like(flux)
            if flux_err is not None:
                return flux.copy(), np.asarray(flux_err, dtype=np.float64).copy(), continuum
            return flux.copy(), None, continuum

    else:
        # ---------------------------------------------------------------
        # Auto-smoothing UnivariateSpline (legacy / fallback)
        # ---------------------------------------------------------------
        w = None
        if flux_err is not None:
            flux_err_arr = np.asarray(flux_err, dtype=np.float64)
            err_sel = flux_err_arr[sel]
            err_sel = np.where(err_sel > 0, err_sel,
                               np.median(err_sel[err_sel > 0])
                               if np.any(err_sel > 0) else 1.0)
            w = 1.0 / err_sel

        spl = UnivariateSpline(wave[sel], flux[sel],
                               k=degree, s=smoothing_factor, w=w)

    continuum = spl(wave)
    continuum = np.where(continuum > 0, continuum, 1.0)

    corrected = flux / continuum

    corrected_err = None
    if flux_err is not None:
        flux_err_arr = np.asarray(flux_err, dtype=np.float64)
        corrected_err = flux_err_arr / continuum

    return corrected, corrected_err, continuum


def _get_plot_wavelength_from_config(config_in):
    """Return the reference wavelength used to choose the diagnostic plot order."""

    spectral_lines = config_in.get('spectral_lines', {})
    if not spectral_lines:
        return None, None

    preferred_label = None
    if 'Na_doublet' in spectral_lines:
        preferred_label = 'Na_doublet'
    else:
        preferred_label = list(spectral_lines.keys())[0]

    lines_dict = spectral_lines.get(preferred_label, {})
    lines_values = list(lines_dict.get('lines', {}).values())

    if len(lines_values) == 0:
        line_range = lines_dict.get('range', None)
        if line_range is None:
            return None, preferred_label
        return 0.5 * (line_range[0] + line_range[1]), preferred_label

    return float(np.mean(lines_values)), preferred_label


def _find_order_for_wavelength(preparation_obs, target_wavelength):
    """Find the order index containing target_wavelength, or the closest order if outside."""

    wave_2d = preparation_obs['wave']
    n_orders = wave_2d.shape[0]

    selected_order = None
    min_distance = np.inf

    for order in range(n_orders):
        wmin = np.nanmin(wave_2d[order, :])
        wmax = np.nanmax(wave_2d[order, :])

        if wmin <= target_wavelength <= wmax:
            return order

        distance = min(abs(target_wavelength - wmin), abs(target_wavelength - wmax))
        if distance < min_distance:
            min_distance = distance
            selected_order = order

    return selected_order


def _median_binned_spline(wave, ratio, bin_size, mask_regions, degree=3,
                         min_bin_points=2):
    """
    Compute a robust broad-band continuum + wiggle model by median-binning
    the spectrum and fitting a cubic spline through the bin medians.

    Bins overlapping ``mask_regions`` are excluded from the spline fit so
    the spline interpolates smoothly through spectral-line regions.

    The median is insensitive to bad pixels / cosmic-ray spikes and the
    bins naturally average over pixel-to-pixel noise.  No polynomial
    pre-detrending is needed.

    Parameters
    ----------
    wave          : 1-D ndarray  – wavelength in Angstrom
    ratio         : 1-D ndarray  – flux ratio (any absolute scale)
    bin_size      : float        – bin width in Angstrom (typical: 2 A)
    mask_regions  : list of [float, float] – wavelength ranges to exclude
    degree        : int          – spline degree (default 3 = cubic)
    min_bin_points: int          – minimum data points per bin to keep it

    Returns
    -------
    continuum : 1-D ndarray or None
        Spline model evaluated at every input wavelength.
        Returns None if the fit cannot be performed.
    """
    finite = np.isfinite(ratio) & (ratio > 0)
    if not np.any(finite):
        return None

    w_min = wave[finite].min()
    w_max = wave[finite].max()
    bin_edges = np.arange(w_min, w_max + bin_size, bin_size)
    if len(bin_edges) < 2:
        return None

    bin_centers = []
    bin_meds    = []
    for i in range(len(bin_edges) - 1):
        b_lo, b_hi = bin_edges[i], bin_edges[i + 1]
        # Exclude bins overlapping any masked region
        if any(b_lo < reg[1] and b_hi > reg[0] for reg in mask_regions):
            continue
        in_bin = finite & (wave >= b_lo) & (wave < b_hi)
        if np.sum(in_bin) < min_bin_points:
            continue
        bin_centers.append(0.5 * (b_lo + b_hi))
        bin_meds.append(float(np.median(ratio[in_bin])))

    if len(bin_centers) < degree + 2:
        return None

    bin_centers = np.array(bin_centers)
    bin_meds    = np.array(bin_meds)

    try:
        spl = UnivariateSpline(bin_centers, bin_meds, k=degree, s=0)
    except Exception:
        return None

    return spl(wave)


def _delete_downstream_caches(config_in, night):
    """
    Remove cached pickle files that were computed from the (uncorrected)
    transmission_preparation so that downstream steps recompute from the
    freshly corrected data.

    Affected patterns per night:
      {output}_{lines}_{night}_transmission_spectrum_*.p
      {output}_{lines}_{night}_transmission_binned_mcmc_*.p
    And global (no-night) aggregates:
      {output}_{lines}_transmission_spectrum_average_*.p
      {output}_{lines}_transmission_binned_mcmc_*.p  (global MCMC results)
    """
    spectral_lines = from_config_get_spectral_lines(config_in)
    output = config_in['output']
    deleted = []

    for lines_label in spectral_lines:
        # Per-night caches
        per_night_patterns = [
            output + '_' + lines_label + '_' + night + '_transmission_spectrum_*.p',
            output + '_' + lines_label + '_' + night + '_transmission_binned_mcmc_*.p',
            # No-lines prefix variant (lines_label may be empty for some configs)
            output + '_' + night + '_transmission_spectrum_*.p',
            output + '_' + night + '_transmission_binned_mcmc_*.p',
        ]
        # Global (no-night) caches
        global_patterns = [
            output + '_' + lines_label + '_transmission_spectrum_average_*.p',
            output + '_' + lines_label + '_transmission_lightcurve_*.p',
            output + '_' + lines_label + '_transmission_binned_mcmc_*.p',
        ]
        for pattern in per_night_patterns + global_patterns:
            for fpath in glob.glob(pattern):
                try:
                    os.remove(fpath)
                    deleted.append(fpath)
                except OSError:
                    pass

    # Note: shared.p is NOT deleted here — it is generated by
    # transmission_spectrum_preparation and does not depend on wiggle correction.

    if deleted:
        print('  Deleted {:d} cached downstream result(s) for night {:s}:'.format(
            len(deleted), night))
        for fpath in deleted:
            print('    - ' + os.path.basename(fpath))
    else:
        print('  No downstream cached results found to delete for night {:s}.'.format(night))


def compute_wiggle_correction(config_in):
    """
    Spectrum-by-spectrum wiggle correction.

    For each observation × echelle order, fit a degree-2 polynomial +
    LSQ spline directly to the individual ``ratio`` array
    (= obs_e2ds / master_out_smooth).

    Because ``ratio`` is already normalised to ~1 in the continuum (the
    master_out_smooth divides away the absolute stellar flux and the Na D
    lines), any residual broad structure and quasi-sinusoidal wiggle in the
    ratio come purely from the individual exposure.  Fitting per-spectrum
    avoids template-mismatch problems caused by temporal variation in wiggle
    amplitude/phase and eliminates the Na D stellar-line contamination that
    affected the master_out approach.

    Algorithm (per obs, per order)
    ------------------------------
    1. Build a boolean mask excluding ``mask_regions`` (Na I D cores) and
       non-finite / non-positive pixels.
    2. Compute ``order_norm = median(ratio[continuum])`` ≈ 1.
    3. Normalise: ``norm_ratio = ratio / order_norm``.
    4. Fit degree-2 polynomial on ``norm_ratio[continuum]`` (broad slope).
    5. Detrend: ``detrended = norm_ratio / poly``.
    6. Fit LSQ B-spline on ``detrended[continuum]`` (wiggle, ~0.5 %).
    7. Apply ``full_model = order_norm × poly × spline`` to
       ratio, ratio_err, deblazed, deblazed_err.

    The Na I D region is excluded from the fit; the spline interpolates
    smoothly through the masked gap, so the planetary absorption signal
    (~0.1 %, much smaller than the ~0.5 % wiggle) is not significantly
    affected.

    Configuration
    -------------
    ::

        wiggle_correction:
          method: spline           # only 'spline' is implemented
          degree: 3                # spline degree (default 3 = cubic)
          knot_spacing: 10.0       # Å between LSQ-spline interior knots
          mask_regions:
            - [5889.0, 5891.0]     # Na I D2
            - [5895.0, 5897.0]     # Na I D1
    """

    night_dict = from_config_get_nights(config_in)

    wiggle_cfg = config_in.get('wiggle_correction', None)
    if wiggle_cfg is None:
        print("{0:45s}  No wiggle_correction section in config – skipping.".format(subroutine_name))
        return

    bin_size     = float(wiggle_cfg.get('bin_size_angstrom', 2.0))
    mask_regions = wiggle_cfg.get('mask_regions', [])

    for night in night_dict:

        try:
            preparation = load_from_cpickle('transmission_preparation', config_in['output'], night)
        except (FileNotFoundError, IOError):
            print("{0:45s} Night:{1:15s}  transmission_preparation not found – skipping.".format(
                subroutine_name, night))
            continue

        if preparation.get('wiggle_corrected_v11', False):
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Already corrected – skipping'))
            continue

        print()
        print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
        print("  bin_size: {0:.1f} A  (median-binned cubic spline)".format(bin_size))
        if mask_regions:
            print("  Masked regions (Angstrom): " + ", ".join(
                "[{:.2f}, {:.2f}]".format(*r) for r in mask_regions))

        lists              = load_from_cpickle('lists',              config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        n_orders = observational_pams['n_orders']
        all_obs  = lists['observations']

        n_corrected = 0
        resid_stds  = []

        for obs in all_obs:
            wave_obs = preparation[obs]['wave']   # (n_orders, n_pixels)

            for order in range(n_orders):
                wave_ord  = wave_obs[order, :]
                ratio_ord = preparation[obs]['ratio'][order, :]

                finite = np.isfinite(ratio_ord) & (ratio_ord > 0)
                if np.sum(finite) < 5:
                    continue

                # Robust median-binned cubic spline continuum model.
                # No polynomial pre-fitting: median binning is inherently
                # robust to bad pixels and captures the broad-band shape
                # (blaze residual, Na D wings) plus the wiggle pattern.
                continuum = _median_binned_spline(
                    wave_ord, ratio_ord,
                    bin_size=bin_size,
                    mask_regions=mask_regions,
                    degree=3,
                    min_bin_points=2)

                if continuum is None:
                    continue

                continuum = np.where(continuum > 0, continuum, 1.0)

                sel = finite & _build_mask(wave_ord, mask_regions)
                if np.any(sel):
                    resid_std = float(np.std((ratio_ord / continuum)[sel]))
                    if resid_std > 1e-6:
                        n_corrected += 1
                        resid_stds.append(resid_std)

                preparation[obs]['ratio'][order, :]        /= continuum
                preparation[obs]['ratio_err'][order, :]    /= continuum
                preparation[obs]['deblazed'][order, :]     /= continuum
                preparation[obs]['deblazed_err'][order, :] /= continuum

        print("  Spectra × orders corrected: {:d} / {:d}".format(
            n_corrected, len(all_obs) * n_orders))
        if resid_stds:
            print("  Residual std after correction (continuum pixels): "
                  "{:.5f}  [{:.5f} – {:.5f}]".format(
                      float(np.mean(resid_stds)),
                      float(np.min(resid_stds)),
                      float(np.max(resid_stds))))

        preparation['wiggle_corrected']     = True
        preparation['wiggle_corrected_v11'] = True
        preparation['wiggle_correction_params'] = {
            'bin_size_angstrom': bin_size,
            'mask_regions':      mask_regions,
        }

        save_to_cpickle('transmission_preparation', preparation, config_in['output'], night)
        _delete_downstream_caches(config_in, night)
        print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Done'))

    print()


def plot_wiggle_correction(config_in, night_input=''):
    """
    Diagnostic plots for the wiggle correction.

    For each night and each in-transit observation, shows:
    - the deblazed transmission spectrum for a representative order  - the fitted spline continuum (wiggles model)
    - the corrected spectrum
    """

    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    wiggle_cfg = config_in.get('wiggle_correction', {})
    mask_regions = wiggle_cfg.get('mask_regions', [])
    target_wavelength, line_label = _get_plot_wavelength_from_config(config_in)

    for night in night_list:

        try:
            preparation = load_from_cpickle('transmission_preparation', config_in['output'], night)
        except (FileNotFoundError, IOError):
            print("No transmission_preparation found for night {}, skipping plots.".format(night))
            continue

        if not preparation.get('wiggle_corrected', False):
            print("Night {}: wiggle correction not applied, no plots available.".format(night))
            continue

        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        # Pick the order containing the target line (Na doublet by default)
        obs_reference = lists['transit_in'][0] if len(lists['transit_in']) > 0 else lists['observations'][0]
        if target_wavelength is not None:
            plot_order = _find_order_for_wavelength(preparation[obs_reference], target_wavelength)
        else:
            plot_order = observational_pams['n_orders'] // 2

        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(lists, observational_pams)

        fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
        if target_wavelength is None:
            fig.suptitle('Wiggle correction - Night: {}  |  Order: {}'.format(night, plot_order))
        else:
            fig.suptitle('Wiggle correction - Night: {}  |  Order: {}  |  {} @ {:.3f} A'.format(
                night, plot_order, line_label, target_wavelength))

        axes[0].set_ylabel('Deblazed (original)')
        axes[1].set_ylabel('Spline continuum')
        axes[2].set_ylabel('Corrected')
        axes[2].set_xlabel('Wavelength [Å]')

        # The stacked continuum is one array per order (not per obs)
        stacked_continuum = preparation.get('wiggle_stacked_continuum', None)
        if stacked_continuum is not None:
            continuum = stacked_continuum[plot_order, :]
        else:
            continuum = np.ones(preparation[obs_reference]['wave'].shape[1])

        axes[1].plot(preparation[obs_reference]['wave'][plot_order, :],
                     continuum, color='black', lw=1.2, label='stacked continuum')

        for obs in lists['transit_in']:
            wave_ord = preparation[obs]['wave'][plot_order, :]
            corrected = preparation[obs]['deblazed'][plot_order, :]

            # Reconstruct original from corrected * continuum (best-effort)
            original = corrected * continuum

            axes[0].plot(wave_ord, original, alpha=0.5, lw=0.8)
            axes[2].plot(wave_ord, corrected, alpha=0.5, lw=0.8)

        # Shade masked regions
        for ax in axes:
            for region in mask_regions:
                ax.axvspan(region[0], region[1], color='gray', alpha=0.2, label='masked')
            ax.axhline(1.0, color='black', lw=0.5, ls='--')

        plt.tight_layout()
        plt.show()
