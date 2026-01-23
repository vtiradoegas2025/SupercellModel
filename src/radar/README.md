# Radar Forward Operators Module

This module implements forward operators that simulate radar observables from atmospheric model state, enabling data assimilation and radar signature validation.

## Overview

Radar forward operators transform prognostic model variables (hydrometeor mixing ratios, wind fields) into radar observables that can be compared with actual radar measurements. This enables:

- **Data assimilation**: Using radar observations to improve model forecasts
- **Operator validation**: Comparing simulated vs. observed radar signatures
- **Model evaluation**: Assessing microphysics and dynamics through radar metrics

## Architecture

```
src/radar/
├── radar.cpp                 # Main radar module coordinator
├── factory.cpp/.hpp          # Scheme factory for radar operators
├── base/
│   └── radar_base.cpp/.hpp   # Common radar utilities and interfaces
└── schemes/                  # Individual radar operator implementations
    ├── reflectivity/         # Z/Ze operators
    ├── velocity/             # Vr operators
    └── zdr/                  # Polarimetric operators
```

## Radar Physics Fundamentals

### Radar Reflectivity (Z/Ze)

Radar reflectivity is proportional to the 6th moment of the drop size distribution (DSD):

```
Z = ∫ N(D) D^6 dD ∝ q_ρ^2 / N_t
```

Where:
- N(D): Drop size distribution (m⁻³ mm⁻¹)
- D: Drop diameter (mm)
- q: Mixing ratio (kg/kg)
- ρ: Hydrometeor density (kg/m³)
- N_t: Total number concentration (m⁻³)

### Equivalent Reflectivity Factor (Z_e)

For Rayleigh scattering (λ >> D), the equivalent reflectivity factor is:

```
Z_e = (λ^4 / π^5 |K|^2) ∫ σ_bscat dV
```

Where:
- λ: Radar wavelength
- K: Dielectric factor (|K|² ≈ 0.93 for water, 0.197 for ice)
- σ_bscat: Backscattering cross-section

### Doppler Radial Velocity (V_r)

The radial component of velocity along the radar beam:

```
V_r = (u_r cosφ + v_θ sinφ) / |k̂| + V_t cosα
```

Where:
- φ: Azimuth angle from radar
- α: Elevation angle
- V_t: Terminal fall velocity
- k̂: Unit vector along radar beam

## Operator Implementations

### Reflectivity Operators

#### Single-Moment Schemes (Kessler, Lin)
```cpp
// Power-law approximation
Z = C * q_r^{1.5-2.0}  // C calibrated for scheme-specific DSD assumptions
```

**Assumptions:**
- Monodisperse or fixed-shape drop size distributions
- Single parameter (mixing ratio) determines reflectivity
- Calibration constants from literature or tuning

#### Double-Moment Schemes (Thompson, Milbrandt)
```cpp
// Explicit DSD integration
Z = ∫ N_0(D) * D^6 * exp(-Λ D) dD
// where N_0 and Λ determined from q_r and N_r
```

**Features:**
- Prognostic number concentrations (N_r, N_i, N_s, N_g)
- More realistic DSD shapes (gamma, exponential)
- Better microphysical consistency

### Velocity Operators

#### Basic Implementation
```cpp
// Radial velocity from winds + fall speeds
V_r[i][j][k] = (u[i][j][k] * cos_phi + v[i][j][k] * sin_phi) / range_factor
             + terminal_fall_speed(q_species[i][j][k], density)
```

**Terminal velocity models:**
- **Rain**: V_t = min(V_tmax, a × q_r^b) with empirical coefficients
- **Ice**: V_t functions of habit, density, and size
- **Hail**: Size-dependent fall speeds with density corrections

### Polarimetric Operators

#### Differential Reflectivity (Z_DR)
```cpp
Z_DR = 10 * log10(Z_H / Z_V)
```

**Physical interpretation:**
- Z_H: Horizontal polarization reflectivity
- Z_V: Vertical polarization reflectivity
- Oblate spheroids (rain) → Z_DR > 0 dB
- Spherical particles → Z_DR ≈ 0 dB

#### Specific Differential Phase (K_DP)
```cpp
K_DP = (180/π) * ∫ (Re(f_HH) - Re(f_VV)) dz / λ
```

**Applications:**
- Liquid water content retrieval
- Rain rate estimation
- Hail detection (negative K_DP)

## Usage Examples

### Basic Reflectivity Calculation
```cpp
#include "radar.hpp"

// Initialize radar operator
auto radar = create_radar_scheme("reflectivity");
radar->initialize(radar_config);

// Compute reflectivity from model state
RadarStateView state{...};  // Model fields
RadarOut output;
radar->compute(radar_config, state, output);

// Access results
float ze_dBZ = output.reflectivity_dbz[i][j][k];
```

### Multi-Operator Pipeline
```cpp
// Initialize multiple operators
auto reflectivity_op = create_radar_scheme("reflectivity");
auto velocity_op = create_radar_scheme("velocity");
auto zdr_op = create_radar_scheme("zdr");

// Apply to same model state
RadarOut ref_out, vel_out, zdr_out;
reflectivity_op->compute(config, state, ref_out);
velocity_op->compute(config, state, vel_out);
zdr_op->compute(config, state, zdr_out);
```

## Configuration

### Radar Operator Settings
```yaml
radar:
  scheme: "reflectivity"  # or "velocity", "zdr"
  operator_tier: "fast_da"  # "fast_da", "psd_moment"
  has_qr: true           # Enable rain reflectivity
  has_qs: true           # Enable snow
  has_qg: true           # Enable graupel
  has_qh: true           # Enable hail
```

### Operator-Specific Parameters
```yaml
radar:
  reflectivity:
    wavelength_mm: 100.0    # S-band radar
    temperature_correction: true
    mixed_phase_handling: "max_reflectivity"
  velocity:
    terminal_velocity_model: "empirical"
    density_correction: true
  zdr:
    canting_angle_model: "gaussian"
    melting_layer_detection: true
```

## Validation and Testing

### Radar Signature Validation
```bash
# Test radar operators against known cases
python tests/validate_simulation.py --radar-validation --timestep 100

# Compare with idealized radar signatures
python tests/radar_validation.py --scheme kessler --case supercell
```

### Operator Consistency Checks
- **Z-V_r relationship**: Expected correlations for different storm types
- **Polarimetric signatures**: Z_DR patterns for various hydrometeors
- **Attenuation effects**: Path-integrated signal loss
- **Beam broadening**: Volume averaging effects

## Implementation Notes

### Performance Considerations
- **Computational cost**: Reflectivity ~ O(N_grid), Velocity ~ O(N_grid)
- **Memory usage**: Additional 3D arrays for radar fields
- **Cache efficiency**: Operators process contiguous data structures

### Numerical Stability
- **Division by zero**: Handle q_r = 0 cases gracefully
- **Range checking**: Validate input field ranges
- **Interpolation**: Smooth transitions at phase boundaries

### Extensibility
- **New operators**: Implement `RadarSchemeBase` interface
- **Hybrid schemes**: Combine multiple operators
- **Advanced physics**: Mie scattering, multiple scattering corrections

## Scientific References

### Reflectivity Operators
- Smith et al. (1975) - Radar reflectivity factor calculations
- Sauvageot (1992) - Rainfall rate relationships
- Thompson et al. (2012) - Volumetric radar sampling operators

### Velocity Operators
- Doviak and Zrnić (1993) - Doppler radar principles
- Miller et al. (1986) - Terminal velocity parameterizations

### Polarimetric Radar
- Zrnić (1987) - Differential reflectivity theory
- Ryzhkov and Zrnić (1998) - Polarimetric applications
- Kumjian (2013) - Polarimetric radar signatures

This module provides the foundation for radar data assimilation and model validation in convective storm research.
