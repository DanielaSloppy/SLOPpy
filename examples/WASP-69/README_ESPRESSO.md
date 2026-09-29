# Using SLOPpy with ESPRESSO data

This guide explains the ESPRESSO-specific behaviour of SLOPpy and how to set up a configuration file for a transit observation.

---

## 1. ESPRESSO data format and DRS products

ESPRESSO is a fiber-fed cross-dispersed echelle spectrograph at the VLT covering 380–788 nm across two CCDs (blue and red).  The ESO Data Reduction Software (DRS) delivers a set of `_S2D_*` FITS files per observation.  The relevant extensions are:

| Extension | Content |
|-----------|---------|
| `[1]` | 2-D flux array (n_orders × n_pixels) |
| `[2]` | 2-D flux error array |
| `[5]` | 2-D wavelength array (barycentric, air) |
| `[7]` | 2-D wavelength step array |

The wavelength axis in the S2D files is given in the **Barycentric rest frame** (BJD_TDB).  SLOPpy shifts the spectra back to the **Observer rest frame** (ORF) on loading, exactly as it does for HARPS/HARPS-N data, so all subsequent pipeline steps work identically.

### 1.1 Odd / even spectral slices

Each spectral order in an ESPRESSO S2D file is stored as **two adjacent rows** corresponding to the two slices of the 2-fiber slit.  SLOPpy treats them as two independent datasets labelled `odd` (rows 1, 3, 5, …) and `even` (rows 0, 2, 4, …).  Define them as separate nights in the YAML file and set:

```yaml
spectral_selection: odd   # or: even
```

Process both slices and average the resulting transmission spectra to maximise S/N.

### 1.2 Available DRS products

| Keyword | File opened | Description |
|---------|-------------|-------------|
| `use_ESO_deblazed: True` (default) | `_S2D_TELL_CORR_A.fits` or `_S2D_SKYSUB_A.fits` or `_S2D_A.fits` | Blaze-divided spectrum; set `calib_data['blaze'] = 1` internally |
| `use_ESO_deblazed: False` | `_S2D_BLAZE_A.fits` | Raw (blazed) spectrum; blaze computed from ratio |
| `use_ESO_telluric_correction: True` | `_S2D_TELL_CORR_A.fits` | Pre-applied telluric + sky correction |
| `use_ESO_sky_correction: True` | `_S2D_SKYSUB_A.fits` | Pre-applied sky correction only |

> **Important – TELL_SPECTRUM wavelength bug**: the `_S2D_TELL_SPECTRUM_A.fits` product has `WAVELENGTH_AIR_BARY` in its header, but the wavelength axis is actually in the **Observer Rest Frame**, not the barycentric frame.  SLOPpy handles this correctly when `apply_ESO_telluric_correction: True`.

### 1.3 Recommended settings

For the vast majority of ESPRESSO analyses use:

```yaml
use_ESO_telluric_correction: True   # read _S2D_TELL_CORR_A.fits
use_ESO_sky_correction:      True   # sky subtraction already applied
use_ESO_deblazed:            True   # no blaze file needed
```

This relies entirely on ESO DRS products and avoids the need for lamp files (which are required for sky correction through SLOPpy but are not yet fully supported).

#### Alternative: apply telluric correction on non-sky-subtracted spectra

```yaml
apply_ESO_telluric_correction: True
use_ESO_sky_correction: False
use_ESO_deblazed: True
```

SLOPpy will read `_S2D_SKYSUB_A.fits` and `_S2D_TELL_CORR_A.fits`, compute the telluric correction spectrum as their ratio, and apply it to `_S2D_A.fits`.  Note that combining `apply_ESO_telluric_correction: True` with `use_ESO_sky_correction: True` is **not allowed** and will cause SLOPpy to abort with an error.

---

## 2. Night definition

```yaml
nights:
  '2022-07-06_odd':
    night: '2022-07-06'        # calendar date of the night (UT)
    instrument: ESPRESSO
    telescope: UT1             # UT1 | UT2 | UT3 | UT4 – required for FITS keywords
    spectral_selection: odd    # odd | even
    use_ESO_telluric_correction: True
    use_ESO_sky_correction: True
    use_ESO_deblazed: True
    all:          2022-07-06_all.list
    in_transit:   2022-07-06_transit_in.list
    out_transit:  2022-07-06_transit_out.list
    telluric_list: 2022-07-06_transit_out.list
    full_transit: 2022-07-06_transit_full.list
    time_of_transit: 2459767.829   # BJD_TDB of mid-transit
    spline_residuals: False
```

The observation list files (e.g. `2022-07-06_all.list`) contain one **filename root** per line (no extension, no path), e.g.:

```
r.ESPRE.2022-07-07T04_47_08.723
r.ESPRE.2022-07-07T04_53_02.029
…
```

SLOPpy appends `_S2D_TELL_CORR_A.fits` (or the relevant suffix depending on the keywords above) and looks for the file in the `data_archive` folder.

---

## 3. Instrument section

```yaml
instruments:
  ESPRESSO:
    data_archive: ../archive/   # path relative to the config file
    orders: red_ccd             # red_ccd | blue_ccd | full_ccd
    mask: K6                    # CCF mask (K6 recommended for K/G-type stars)
    wavelength_rescaling: [5869.0, 5877.0]
    resolution: 140000          # HR mode; use 70000 for MR mode
```

`orders` selects which CCD orders to load:

| Value | Orders | Wavelength range (approx.) |
|-------|--------|---------------------------|
| `red_ccd` | 47–84 | ~5200–7900 Å |
| `blue_ccd` | 0–46 | ~3800–5200 Å |
| `full_ccd` | 0–84 | full range |

---

## 4. Wiggle correction (ESPRESSO-specific)

ESPRESSO transmission spectra show quasi-sinusoidal oscillations (wiggles) in the continuum, believed to originate from thin-film interference in the fiber train and dichroic.  These must be corrected **after** computing the transmission spectrum ratio but **before** the CLV+RM correction and the MCMC fit.

SLOPpy implements a cubic-spline correction following the approach of Yan et al. (2021, A&A 645, A22):

1. For each observation and each spectral order, a cubic `UnivariateSpline` (k = 3) is fitted to the continuum of the deblazed transmission ratio, **excluding the spectral-line regions** listed in `mask_regions`.
2. The spectrum is divided by the fitted spline, removing the wiggles without distorting the planetary absorption lines.

### 4.1 Pipeline placement

```yaml
pipeline:
  - differential_refraction
  - master_out
  - transmission_spectrum_preparation
  - wiggle_correction             # ← here, after prep, before CLV+RM / MCMC
  - transmission_binned_mcmc
  - transmission_spectrum
  - transmission_spectrum_average
```

### 4.2 Configuration

```yaml
wiggle_correction:
  method: spline          # only 'spline' is currently implemented
  degree: 3               # cubic spline (recommended)
  smoothing_factor: null  # null → scipy automatic selection (recommended)
                          # increase (e.g. 0.05–0.5) for stronger smoothing
  mask_regions:           # wavelength intervals [Å] to exclude from the fit
    - [5889.0, 5891.0]    # Na I D2
    - [5895.0, 5897.0]    # Na I D1
    # Add more regions if other lines are present in the spectral window
```

`smoothing_factor` is the `s` parameter of `scipy.interpolate.UnivariateSpline`.  Setting it to `null` lets scipy select `s = N` (number of data points), which works well for typical ESPRESSO data.  If the correction is too aggressive (removing real signal) increase `s`; if it leaves residuals, decrease it.

### 4.3 Diagnostic plot

Add `wiggle_correction` to the `plots` section of the config to display, for each in-transit observation, a three-panel figure showing:
- the original deblazed transmission ratio per order
- the fitted spline continuum (the wiggle model)
- the corrected spectrum

### 4.4 Method comparison

| Method | Pros | Cons |
|--------|------|------|
| **Cubic spline (iSpec, this implementation)** | Masks line regions; flexible smoothing; per-order 2-D correction | Smoothing factor must be tuned |
| Sinusoidal fit (ANTARESS) | Physically motivated for interference patterns | Risk of fitting planetary signal at the same frequency; harder to generalise |
| Bivariate spline | Works on 2-D (wave × time) data | Overkill for 1-D per-order spectra; harder to control |

The cubic-spline approach is preferred because it makes no assumption about the wiggle frequency and safely avoids the planetary absorption features.

---

## 5. Sky correction

Sky correction through SLOPpy (`sky_correction` step) requires lamp files which are not currently available from the standard ESO DRS ESPRESSO archive.  **Use the ESO pre-subtracted sky products instead** (`use_ESO_sky_correction: True`).

---

## 6. Complete pipeline for ESPRESSO

A recommended minimal pipeline for ESPRESSO transit spectroscopy is:

```yaml
pipeline:
  - differential_refraction
  - master_out
  - transmission_spectrum_preparation
  - wiggle_correction
  - clv_rm_models_lines
  - transmission_binned_mcmc
  - transmission_spectrum
  - transmission_spectrum_average
```

Optional steps:
- `telluric_molecfit_preparation` + `telluric_molecfit_coadd` – if a dedicated Molecfit telluric correction is needed beyond the ESO DRS product.
- `telluric_template` – template-based telluric correction.
- `differential_refraction_update` – second-pass differential-refraction correction.

---

## 7. Example files

The `examples/WASP-69/` directory contains:

| File | Description |
|------|-------------|
| `WASP-69_transit_ESPRESSO.yaml` | Complete configuration file |
| `2022-07-06_all.list` | Full list of observations for the night 2022-07-06 |
| `2022-07-06_transit_in.list` | In-transit observations *(to be created by the user)* |
| `2022-07-06_transit_out.list` | Out-of-transit observations *(to be created by the user)* |
| `2022-07-06_transit_full.list` | Full-transit observations *(to be created by the user)* |

Raw ESPRESSO data for WASP-69b (night 2022-07-06) can be downloaded from:  
https://drive.google.com/drive/folders/1wlgUsp-EHsHEQ3B1lL5tsmhP0zMlhPb8

Place the downloaded `_S2D_TELL_CORR_A.fits` files (and any other required products) in a folder named `archive/` at the same level as the `examples/WASP-69/` directory (or update `data_archive` accordingly).

---

## 8. References

- Allart et al. 2017, A&A 606, A144 – first ESPRESSO Na detection in HD 189733b
- Casasayas-Barris et al. 2021, A&A 647, A26 – ESPRESSO transmission spectroscopy
- Yan et al. 2021, A&A 645, A22 – iSpec-based cubic spline continuum correction
- Bourrier et al. (ANTARESS) – sinusoidal wiggle correction approach
- ESO ESPRESSO DRS User Manual – description of S2D product formats
