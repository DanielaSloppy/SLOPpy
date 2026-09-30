# Emission Spectra (Eclipse Nights)

SLOPpy was originally designed to extract **transmission spectra**, observed while the
planet transits in front of its host star. It can now also extract **emission
spectra**, observed when the planet does **not** transit the star (e.g. close to the
secondary eclipse, when the day side is visible).

A `phase` keyword marks each night as a transit or an eclipse night, and a single
configuration file can contain nights of both types. Emission spectra have their own
steps, parallel to the transmission ones:

| Transmission (`phase: transit`)     | Emission (`phase: eclipse`)     |
|-------------------------------------|---------------------------------|
| `transmission_spectrum_preparation` | `emission_spectrum_preparation` |
| `transmission_spectrum`             | `emission_spectrum`             |
| `transmission_spectrum_average`     | `emission_spectrum_average`     |
| `write_output_transmission`         | `write_output_emission`         |

Each step analyzes only the nights of its own type and writes its own output files,
so transit and eclipse nights are never mixed when spectra from several nights are
averaged.

## Configuration

### The `phase` keyword

Add the `phase` keyword to each night in the `nights` section of the YAML file, and
list the steps you need in the `pipeline` section:

```yaml
pipeline:
  - ...
  - master_out
  - transmission_spectrum_preparation   # nights with phase: transit
  - emission_spectrum_preparation       # nights with phase: eclipse
  - emission_spectrum                   # one per spectral line set, each eclipse night
  - emission_spectrum_average           # one per spectral line set, all eclipse nights
  - write_output_emission               # full spectrum, each eclipse night
plots:
  - emission_spectrum_preparation
  - emission_spectrum
  - emission_spectrum_average
  - write_output_emission
nights:
  '2017-10-24':
    all: 2017-10-24_all.list
    in_transit: 2017-10-24_transit_in.list
    out_transit: 2017-10-24_transit_out.list
    full_transit: 2017-10-24_transit_full.list
    instrument: HARPS
    mask: G2
    time_of_transit: 2458051.6680533197
    phase: eclipse
  '2017-11-22':
    ...
    phase: transit
```

If `phase` is missing, it defaults to `transit`, so existing configuration files
keep working without changes. The values are case-sensitive: any value other than
`transit` or `eclipse` stops the pipeline with a `ValueError`.

Steps that run before the preparation (sky correction, telluric correction,
differential refraction, master-out, and so on) still process every night,
whatever its `phase`.

### Time of the secondary eclipse

For each eclipse night, SLOPpy needs the central time of the secondary eclipse,
`Tc,ecl`. You can give it in either of two ways:

| Night keywords                        | How `Tc,ecl` is obtained |
|---------------------------------------|--------------------------|
| `time_of_transit` only                | Computed from the transit: `time_of_transit + P/2` for circular orbits, or the Keplerian solution for eccentric orbits (see below). Among the eclipses separated by one orbital period, SLOPpy picks the one closest to the observations of the night, so the transit can be any epoch, before or after the night |
| `time_of_eclipse`                     | Used as given. If `time_of_transit` is missing, it is derived from `time_of_eclipse` by moving back to the preceding transit |

```yaml
nights:
  '2017-10-24':
    ...
    phase: eclipse
    time_of_transit: 2458051.6680533197   # Tc,ecl = time_of_transit + P/2 (circular orbit)
  '2017-10-26':
    ...
    phase: eclipse
    time_of_eclipse: 2458054.3829         # Tc,ecl given explicitly
```

`time_of_transit` always keeps its meaning of transit time. The planet radial
velocities are computed from it, so it must refer to the transit even for eclipse
nights.

For **eccentric orbits** (`orbit: eccentric`), the eclipse does not fall exactly
half an orbit after the transit. The time between transit and eclipse is computed
from the mean anomalies of the two conjunctions. The true anomaly is `π/2 - ω` at the
transit and `3π/2 - ω` at the eclipse, with the same ω convention that SLOPpy uses
for the transit time and the planet radial velocity. The orbital inclination is
neglected in this computation.

### Observation lists for eclipse nights

Eclipse nights use the same list keywords as transit nights, with the following
roles:

| Keyword        | Role for eclipse nights | Automatic selection |
|----------------|-------------------------|---------------------|
| `out_transit`  | Reference spectra, used to build the master-out | Planet completely behind the star: the whole exposure falls between the second and third contact of the eclipse |
| `full_transit` | Spectra averaged to obtain the emission spectrum | Planet completely outside the stellar disk: the whole exposure falls before the first or after the fourth contact |
| `in_transit`   | Spectra shown in the plots of `emission_spectrum_preparation` | All observations not in `out_transit`, i.e. with the planet at least partially in view (including ingress and egress) |

As for transit nights, the lists are built automatically when the list files do not
exist, and they are then written to the files named in the configuration. If the
files already exist, for example from a previous analysis of the same night as a
transit night, they are used as they are. Delete them to build the lists again.

The whole exposure (`BJD ± EXPTIME/2`) is checked, so the master-out never
includes a spectrum taken while the planet was partially visible. If no
observation falls entirely within the full eclipse, SLOPpy prints a warning,
because the master-out cannot be computed.

The windows use the durations of the eclipse, taken from the `planet` section:

| Keyword                  | Meaning                                         | If missing |
|--------------------------|-------------------------------------------------|------------|
| `total_eclipse_duration` | First to fourth contact of the eclipse [days]   | `total_transit_duration`, rescaled for eccentric orbits |
| `full_eclipse_duration`  | Second to third contact of the eclipse [days]   | `full_transit_duration`, rescaled for eccentric orbits |

For circular orbits, eclipse and transit have the same durations. For eccentric
orbits, the transit durations are multiplied by `(1 + e sin ω) / (1 - e sin ω)`
(Winn 2010, Eq. 16). If only `transit_duration` is available, it is used for both
windows. SLOPpy then prints a warning, because ingress and egress observations
could enter the master-out.

### Stellar radial velocity

Unless `use_rv_from_ccf` or `use_analytical_rvs` is set, the stellar radial
velocity is modelled as a straight line fitted to the RVs of the `out_transit`
spectra. The intercept is used as the systemic velocity. For transit nights the
line is centred on `time_of_transit`, as before. For eclipse nights it is centred
on `Tc,ecl`: centring it on the transit would extrapolate the fit by half an orbit
and bias the systemic velocity. The time used is stored in
`observational_pams['time_of_RV_reference']`. The eclipse time is stored in
`observational_pams['time_of_eclipse']`.

### Composite master-out

With `use_composite: True` in the `master-out` section, the composite master-out is
computed separately for transit and eclipse nights. Eclipse nights use a composite
built only from their own `out_transit` spectra, taken with the planet behind the
star. Transit nights use a composite built only from the other transit nights. The
composite stores its `phase`. Composites saved by earlier versions of SLOPpy have
no `phase` and are treated as transit composites.

### Reference frames

`emission_spectrum`, `emission_spectrum_average` and `write_output_emission` compute
the spectrum in the planet rest frame. The `_planetRF`, `_stellarRF` and
`_observerRF` versions select the frame explicitly (e.g.
`emission_spectrum_stellarRF`). The planet-frame shift uses the planet radial
velocity at the time of each observation, so it applies to any orbital phase.

### Continuum normalization parameters

In `emission_spectrum_preparation`, both the observed spectra and the master-out are
divided by their continuum before subtraction. The continuum is estimated order by
order with a Savitzky-Golay filter. Three optional keywords control the filter:

| Keyword             | Default | Description                                                              |
|---------------------|---------|--------------------------------------------------------------------------|
| `savgol_window`     | `301`   | Width of the filter window, in pixels                                    |
| `savgol_polyorder`  | `2`     | Order of the polynomial fitted within the window                         |
| `savgol_window_kms` | none    | Width of the filter window, in km/s. When set, it overrides `savgol_window` |

A window defined in pixels covers different velocity ranges on different
instruments, because each instrument samples a resolution element with a
different number of pixels. For example, 301 pixels cover a wider velocity range
on HARPS than on ESPRESSO. `savgol_window_kms` gives the same physical width on
every instrument: for each order, SLOPpy converts it into pixels using the local
velocity step of that order:

```
window_pixels = window_kms / (c * median(step / wave))
```

Whichever keyword you use, SLOPpy adjusts the window if needed: it is made odd,
kept larger than `savgol_polyorder`, and kept no longer than the order.

You can set these keywords for a single night or in the `instruments` section. A
value in the `instruments` section applies to every night observed with that
instrument, unless the night sets its own value:

```yaml
nights:
  '2017-10-24':
    ...
    phase: eclipse
    savgol_window_kms: 250.   # this night only
instruments:
  ESPRESSO:
    ...
    savgol_window: 501        # every ESPRESSO night without its own value
    savgol_polyorder: 2
```

## Method

### Preparation

For every observation, the master-out is first moved to the observer reference
frame and rebinned onto the 2D (e2ds) wavelength grid of the observation. This
step is the same for transit and eclipse nights.

**Transmission** (`transmission_spectrum_preparation`, unchanged):

```
ratio     = F / M
ratio_err = ratio * sqrt( (σ_F / F)^2 + (σ_M / M)^2 )
```

where `F` is the e2ds spectrum of the observation and `M` the rebinned master-out.
The result is then divided by the blaze function and by the normalized pixel step.

**Emission** (`emission_spectrum_preparation`):

```
F_norm    = F / C_F
M_norm    = M / C_M
ratio     = F_norm - M_norm
ratio_err = sqrt( (σ_F / C_F)^2 + (σ_M / C_M)^2 )
```

where `C_F` and `C_M` are the Savitzky-Golay continua of the observation and of
the master-out. Continuum values equal to zero are replaced by one. The
uncertainty of the continuum itself is not propagated.

The results are stored under the same keys used for transmission spectra
(`ratio`, `ratio_err`, `deblazed`, `deblazed_err`, `ratio_precleaning`, ...).
Unlike transmission spectra, emission spectra are centred on **zero** rather than
on unity.

Other differences from the transmission preparation:

- **No deblazing.** The continuum normalization already removes the blaze
  function and the pixel-step dependence, so `deblazed` is a copy of `ratio`.
- **Telluric residual cleaning** (`spline_residuals: True`). For transmission, each
  spectrum is divided by its median before the spline fit. That would mean
  dividing by a value close to zero, so for emission the spectra are shifted by +1
  instead, and the shift is removed afterwards. The 5% rejection threshold then
  corresponds to an absolute threshold of 0.05 on the emission spectrum.
- **No rescaling to unity.** The `rescaling` and `rescaled` keys are not computed,
  even with `full_output: True`.
- **Plots.** The emission maps are shifted by +1, so that they can share the colour
  scale of the transmission maps. The spectra of the `in_transit` observations are
  plotted without any shift.

### Emission spectrum of each night

`emission_spectrum` works on the wavelength range of each set of spectral lines
(the `range` keyword of the `spectral_lines` section). `write_output_emission` does
the same over the full spectrum, or over the `range` of the `full_spectrum`
section when one is given.

For each observation, both steps:

1. shift the 2D emission spectrum to the selected reference frame;
2. rebin it to the 1D wavelength grid shared by all nights;
3. average the spectra of the `full_transit` list, weighting each pixel by the
   inverse of its variance. They do the same for the `out_transit` list;
4. rebin both averages to the coarser binned grid.

Pixels not covered by any observation are excluded from the average. In the output,
they have both value and error equal to zero.

`transmission_spectrum` and `write_output_transmission` apply three more
corrections. None of them applies to emission spectra:

| Correction in the transmission steps  | Why it is not applied to emission spectra |
|---------------------------------------|-------------------------------------------|
| Division by the blaze function during the rebinning | The blaze has already been removed by the continuum normalization |
| CLV+RM correction (division by the model) | There is no Rossiter-McLaughlin effect or centre-to-limb variation outside the transit |
| Continuum normalization (division by a polynomial or Savitzky-Golay fit) | The spectra are already continuum-normalized, and a division is meaningless for spectra centred on zero |

The average of the `out_transit` list works as a null test. It should be
consistent with zero: a residual signal there reveals problems in the master-out
or in the continuum normalization, for example near the edges of the orders.

Unlike `transmission_spectrum`, no MCMC results are used. The output files carry
no solution suffix: `emission_spectrum_planetRF` rather than
`transmission_spectrum_planetRF_user`.

### Average of the eclipse nights

`emission_spectrum_average` loads the emission spectra computed by
`emission_spectrum` for every eclipse night. It then averages all the individual
observations of all the nights, again weighting each pixel by the inverse of its
variance. It does this separately for the `full_transit` and `out_transit` lists.
If the spectrum of any night is missing, the average is skipped.

The average is always recomputed, as it only combines already computed spectra.
The spectra of each night are also stored in the output, under the name of the
night, and plotted below the global average.

### Output files

| Step                             | Output file (`<output>_[<lines>_][<night>_]<name>.p`) |
|----------------------------------|-------------------------------------------------------|
| `emission_spectrum_preparation`  | `emission_preparation`, one per night                 |
| `emission_spectrum`              | `emission_spectrum_<reference>`, one per night and set of lines |
| `emission_spectrum_average`      | `emission_spectrum_average_<reference>`, one per set of lines |
| `write_output_emission`          | `write_output_emission_<reference>`, one per night    |

Main keys of `emission_spectrum` and `write_output_emission` outputs:

| Key                                        | Content                                                         |
|--------------------------------------------|-----------------------------------------------------------------|
| `wave`, `step`                             | 1D wavelength grid                                              |
| `<obs>['rebinned']`, `<obs>['rebinned_err']` | Emission spectrum of each observation, in the selected frame  |
| `average`, `average_err`                   | Weighted average of the `full_transit` spectra                  |
| `binned`, `binned_err`                     | `average` rebinned to `binned_wave`                             |
| `average_out`, `binned_out` (and `_err`)   | Same, for the `out_transit` spectra                             |
| `savgol_pams`                              | Normalization parameters used by the preparation                |

`emission_preparation` and `emission_spectrum` store `savgol_pams`, the
normalization parameters used by the preparation. When you rerun a step, a stored
result is reused only if these parameters still match. If you change any
`savgol_*` keyword, both steps are computed again automatically.

### Plots

Plots are saved in the `plots` folder:

- `emission_spectrum` and `write_output_emission`, for each night: the spectra of
  the individual observations (`_observations.pdf`), and the averages with their
  binned versions (`_binned.pdf`). Each figure has one panel for the `full_transit`
  list and one for the `out_transit` list.
- `emission_spectrum_average`: the average over all nights, with the average of
  each night plotted below it (`_binned.pdf`).

The vertical limits are symmetric around zero and adapt to the data.

## Transit-only steps

The steps that read `transmission_preparation` now process only the nights with
`phase: transit` (or no `phase` keyword). Their plot functions skip eclipse nights
too. These steps are:

- `transmission_spectrum`, `transmission_spectrum_average`
- `transmission_binned_mcmc`, `transmission_mcmc`
- `write_output_transmission`
- `wiggle_correction`, `sysrem_correction`
- `quick_transmission`

Other transit-specific steps, such as `transmission_lightcurve`,
`spectra_lightcurve` and their `_average` versions, still process every night. Do
not use them on eclipse nights.

## Code changes

| File | Change |
|------|--------|
| `SLOPpy/emission_spectrum_preparation.py` | New step: `compute_emission_spectrum_preparation()`, `plot_emission_spectrum_preparation()`, and the continuum function `compute_savgol_continuum()` |
| `SLOPpy/emission_spectrum.py` | New step: `compute_emission_spectrum()` and `plot_emission_spectrum()`, with `_planetRF`, `_stellarRF` and `_observerRF` versions. The computation is shared with `write_output_emission`: see `compute_emission_spectrum_core()`, `plot_emission_spectrum_core()` and `weighted_average_spectra()` |
| `SLOPpy/emission_spectrum_average.py` | New step: `compute_emission_spectrum_average()` and `plot_emission_spectrum_average()`, with the reference-frame versions |
| `SLOPpy/write_output_emission.py` | New step: `write_output_emission()` and `plot_output_emission()`, with the reference-frame versions |
| `SLOPpy/subroutines/preparation_subroutines.py` | New module with the code shared by the two preparation steps: `load_master_out_for_preparation()` (composite and smoothed master-out options) and `rebin_master_out_to_observation()` (shift to the ORF, rebinning to the 2D grid, replacement of bad values) |
| `SLOPpy/transmission_spectrum_preparation.py` | Uses the shared subroutines and skips eclipse nights. The computation of transmission spectra is unchanged |
| `SLOPpy/subroutines/kepler_exo.py` | New functions `kepler_eclipse_time_offset()` (time between transit and secondary eclipse, circular and eccentric orbits) and `kepler_eclipse_duration_ratio()` (ratio of eclipse to transit duration, Winn 2010) |
| `SLOPpy/prepare_datasets.py` | Automatic lists for eclipse nights (`_write_eclipse_list()`, `_get_time_of_eclipse()`, `_get_eclipse_durations()`). For eclipse nights, the linear fit of the stellar RVs is centred on the eclipse. The phase of each night is printed |
| `SLOPpy/master_out.py` | Separate composite master-outs for transit and eclipse nights |
| `SLOPpy/subroutines/io_subroutines.py` | New function `get_night_phase(night_dict, night)`. `from_config_get_nights()` takes an optional `phase` argument that selects the nights of that phase. Without it, all nights are returned, as before. For eclipse nights with `time_of_eclipse` and no `time_of_transit`, the configuration parser derives `time_of_transit` |
| The transit-only steps listed above | Call `from_config_get_nights(config_in, phase='transit')`. Their plot functions skip eclipse nights |
| `SLOPpy/config_default.py` | `savgol_window`, `savgol_polyorder` and `savgol_window_kms` added to `copy_from_instrument`, so they can be set in the `instruments` section |
| `SLOPpy/__init__.py`, `SLOPpy/sloppy_run.py` | Registration of the new steps in the `pipeline` and `plots` sections |

## Possible developments

- Light curves of emission lines, analogous to `transmission_lightcurve`.
- Fit of the emission lines, analogous to `transmission_binned_mcmc`.
- Optional removal of a residual baseline from the emission spectra, by
  subtracting a low-order polynomial.
