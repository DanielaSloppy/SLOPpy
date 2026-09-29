# Eccentric Orbit Support in SLOPpy

Starting from version 1.5.1, SLOPpy fully supports eccentric planetary orbits 
for computing the planet's radial velocity and the CLV+RM correction.

## Configuration

### Circular Orbit (Default)

For circular orbits (eccentricity ≈ 0), use the standard configuration:

```yaml
planet:
  orbit: circular
  period: [2.21857567, 0.000015]                    # Orbital period [days]
  RV_semiamplitude: [151.35, 10.0]                  # Planet RV semi-amplitude [km/s]
  transit_duration: [0.0760, 0.0017]                # Transit duration [days]
  inclination: [85.71, 0.02]                        # Orbital inclination [degrees]
  reference_time_of_transit: [2454279.436714, 0.000015]  # Transit time [BJD]
  radius_ratio: [0.15444, 0.00519]                  # Radius ratio Rp/Rs
  semimajor_axis_ratio: [8.84, 0.27]                # Semi-major axis a/Rs
  impact_parameter: [0.6631, 0.0023]                # Impact parameter
```

### Eccentric Orbit

For orbits with significant eccentricity (e > 0.01), add the following parameters:

```yaml
planet:
  orbit: eccentric                                   # IMPORTANT: specify 'eccentric'
  period: [4.2345678, 0.000012]                      # Orbital period [days]
  RV_semiamplitude: [180.5, 15.0]                    # Planet RV semi-amplitude [km/s]
  eccentricity: [0.25, 0.02]                         # Orbital eccentricity
  omega: [85.3, 2.5]                                 # Argument of periastron [degrees]
  transit_duration: [0.0890, 0.0020]                 # Transit duration [days]
  inclination: [87.2, 0.3]                           # Orbital inclination [degrees]
  reference_time_of_transit: [2459123.456789, 0.00002]  # Transit time [BJD]
  reference_time: 2459120.0                          # Reference time for ephemeris [BJD]
  radius_ratio: [0.12345, 0.00321]                   # Radius ratio Rp/Rs
  semimajor_axis_ratio: [12.5, 0.5]                  # Semi-major axis a/Rs
  impact_parameter: [0.45, 0.05]                     # Impact parameter
```

## Orbital Parameters

### Common Parameters (Circular and Eccentric)

| Parameter | Description | Unit | Format |
|-----------|-------------|------|--------|
| `orbit` | Orbit type | - | `'circular'` or `'eccentric'` |
| `period` | Orbital period | days | [value, error] |
| `RV_semiamplitude` | Planet RV semi-amplitude | km/s | [value, error] |
| `inclination` | Orbital inclination | degrees | [value, error] |
| `reference_time_of_transit` | Transit center time | BJD | [value, error] |
| `radius_ratio` | Radius ratio Rp/Rs | - | [value, error] |
| `semimajor_axis_ratio` | Semi-major axis a/Rs | - | [value, error] |

### Additional Parameters for Eccentric Orbits

| Parameter | Description | Unit | Format | Default |
|-----------|-------------|------|--------|---------|
| `eccentricity` | Orbital eccentricity | - | [value, error] | Required if orbit='eccentric' |
| `omega` | Argument of periastron | degrees | [value, error] | Required if orbit='eccentric' |
| `reference_time` | Reference epoch $T_{ref}$ for the Keplerian time axis | BJD | scalar | **Optional**: if omitted, defaults to `reference_time_of_transit[0]` |

> **Note on time keywords:**
> - `time_of_transit` in the **night** section is the transit center time for that specific observed night.
>   It is used by `prepare_datasets` to classify in/out-of-transit observations and compute RV shifts.
> - `reference_time_of_transit` in the **planet** section is the ephemeris reference epoch $T_0$.
>   It is used by the CLV+RM modules to compute the transit offset along the Keplerian time axis:
>   `Tcent_Tref = reference_time_of_transit[0] - reference_time`.
> - `reference_time` in the **planet** section sets the zero point $T_{ref}$ of the Keplerian time axis.
>   Observation times are shifted as `BJD - T_{ref}` before being passed to Keplerian functions.

## Technical Notes

### Planet Radial Velocity Computation

For circular orbits (e < 10⁻³), the radial velocity is computed using the formula:

```
RV_planet = -K_planet * sin(2π * (t - Tc) / P)
```

For eccentric orbits, the full Keplerian solution is used, which solves 
Kepler's equation for the eccentric anomaly and then computes the true anomaly.

### CLV+RM Correction

The CLV+RM correction automatically accounts for eccentricity when computing 
the planet's position on the stellar disk during transit. For eccentric 
orbits, the orbital distance varies according to:

```
r(f) = a * (1 - e²) / (1 + e * cos(f))
```

where `f` is the true anomaly.

### Backward Compatibility

The code is fully backward compatible. If the `orbit` parameter is not 
specified or is set to `'circular'`, the behavior is identical to 
previous versions.

### MCMC Support

The eccentric orbit support extends to the MCMC transmission spectrum fitting
routines (`transmission_mcmc` and `transmission_binned_mcmc`). When the MCMC
sampler derives updated values for the planet's RV semi-amplitude (K_planet),
the RV shift calculations use the eccentric orbit formulas if specified in
the configuration. The orbital parameters (eccentricity, omega) are read from
the `observational_pams['orbital_parameters']` dictionary set during the
`prepare_datasets` step.

## Usage Examples

### Example 1: WASP-14b (e ≈ 0.09)

```yaml
planet:
  orbit: eccentric
  period: [2.2437661, 0.0000011]
  RV_semiamplitude: [991.0, 15.0]
  eccentricity: [0.091, 0.003]
  omega: [252.3, 1.5]
  inclination: [84.32, 0.10]
  reference_time_of_transit: [2454463.57583, 0.00020]
  radius_ratio: [0.0991, 0.0007]
  semimajor_axis_ratio: [5.711, 0.025]
```

### Example 2: HD 80606b (e ≈ 0.93)

```yaml
planet:
  orbit: eccentric
  period: [111.4367, 0.0004]
  RV_semiamplitude: [465.0, 20.0]
  eccentricity: [0.9330, 0.0012]
  omega: [300.77, 0.10]
  inclination: [89.285, 0.023]
  reference_time_of_transit: [2455210.6430, 0.0016]
  reference_time: 2455210.0  # Use a time close to transit for numerical stability
  radius_ratio: [0.1001, 0.0006]
  semimajor_axis_ratio: [96.0, 2.0]
```

## Troubleshooting

### Error: KeyError 'eccentricity'

If you specify `orbit: eccentric` but `eccentricity` and/or `omega` are missing, 
an error will be raised. Make sure all required parameters are present.

### Unexpected RV Values

For highly eccentric orbits (e > 0.8), verify that:
1. `omega` is specified in degrees (not radians)
2. `reference_time` is close to the observation time to avoid error accumulation
3. All times are in the same system (BJD_TDB recommended)

## References

- Murray & Correia (2010) - "Keplerian Orbits and Dynamics of Exoplanets"
- Sicilia et al. (2022) - Original SLOPpy paper

---

*Documentation updated: April 2026*
