# Sub-Grid Turbulence Module

This module implements parameterizations for unresolved turbulent motions in large-eddy simulations (LES) and convection-permitting models.

## Overview

Sub-grid scale (SGS) turbulence parameterizations account for the effects of eddies smaller than the model grid spacing. These schemes provide:

- **Energy dissipation**: Conversion of kinetic energy to heat
- **Momentum mixing**: Vertical transport of horizontal momentum
- **Scalar transport**: Mixing of heat, moisture, and other tracers
- **Length scales**: Turbulent mixing based on grid resolution

## Architecture

```
src/turbulence/
├── turbulence.cpp            # Module coordinator
├── factory.cpp/.hpp          # Scheme factory
├── base/
│   └── eddy_viscosity.cpp/.hpp # Common turbulence utilities
└── schemes/                  # SGS implementations
    ├── smagorinsky/          # Eddy-viscosity model
    └── tke/                  # Prognostic TKE
```

## Turbulence Fundamentals

### Eddy Viscosity Concept

Turbulent stresses are parameterized using an eddy viscosity ν_t:

```
-ρ <u_i' u_j'> = ρ ν_t (∂U_i/∂x_j + ∂U_j/∂x_i - (2/3)δ_ij ∂U_k/∂x_k)
```

Where ν_t is computed from SGS kinetic energy and length scales.

### Turbulent Kinetic Energy (TKE)

The prognostic TKE equation describes SGS energy evolution:

```
d(tke)/dt = P + B - ε + T + transport
```

Where:
- **P**: Shear production (-u_i'u_j' ∂U_i/∂x_j)
- **B**: Buoyancy production (β θ' w')
- **ε**: Dissipation (c_ε tke^{3/2}/l)
- **T**: Transport terms
- **β**: Buoyancy parameter (g/θ_0)

## Scheme Implementations

### Smagorinsky-Lilly Eddy Viscosity

**Concept**: Local equilibrium between production and dissipation

**Eddy viscosity:**
```cpp
ν_t = (C_s Δ)² √(2 S_ij S_ij) × Φ(z/L)
```

Where:
- C_s = 0.18 (Smagorinsky constant)
- Δ = (Δx Δy Δz)^{1/3} (grid length scale)
- S_ij = (1/2)(∂U_i/∂x_j + ∂U_j/∂x_i) (strain rate tensor)
- Φ(z/L) = stability correction function

**Features:**
- **Algebraic**: No prognostic equations
- **Local**: Depends only on local gradients
- **Stable**: Conservative and numerically robust

### Prognostic TKE Scheme

**TKE Equation:**
```cpp
∂(e)/∂t = -u_i'∂U_i/∂x_j u_j' + β θ' w' - c_ε e^{3/2}/l + ∂/∂x_j [(ν_t + ν) ∂e/∂x_j]
```

**Length Scale Determination:**
```cpp
l = min(l_max, κ z / (1 + κ z / l_max))  // Near-wall
l = min(l_max, c_l √(e)/N)               // Buoyancy-limited
```

Where:
- e: SGS TKE (m²/s²)
- l: Mixing length (m)
- N: Brunt-Väisälä frequency (s⁻¹)
- κ = 0.4: von Kármán constant

**Eddy Viscosity:**
```cpp
ν_t = c_μ √(e) l
```

With c_μ = 0.1 calibrated for neutral conditions.

## Turbulence Closure

### Near-Wall Treatment

**Law of the wall** for surface layer:
```cpp
u(z) = (u_*/κ) ln(z/z0)    // Logarithmic profile
```

**Mixing length** in surface layer:
```cpp
l = κ z (1 - z/δ)²         // Damped near surface
```

### Stability Corrections

**Monin-Obukhov similarity** for stable/unstable regimes:
```cpp
Φ_m(ζ) = 1 + β ζ for stable (ζ > 0)
Φ_m(ζ) = (1 - γ ζ)^{-1/4} for unstable (ζ < 0)
```

Where ζ = z/L (stability parameter), L = -u_*^3 θ_0 / (κ g Q_0)

## Configuration

### Scheme Selection
```yaml
turbulence:
  scheme: "smagorinsky"    # smagorinsky, tke
  dt_sgs: 1.0              # SGS time step (s)
  Cs: 0.18                 # Smagorinsky constant
  Pr_t: 1.0                # Turbulent Prandtl number
```

### Advanced Parameters
```yaml
turbulence:
  smagorinsky:
    stability_correction: true
    dynamic_procedure: false    # Dynamic Smagorinsky-Lilly

  tke:
    c_mu: 0.1                  # Eddy viscosity constant
    c_eps: 0.7                 # Dissipation constant
    l_max: 100.0               # Maximum mixing length (m)
    e_min: 1e-6                # Minimum TKE (m²/s²)
```

## Implementation Details

### Grid-Level Application

**Momentum tendencies:**
```cpp
du/dt = ∂/∂z (ν_t ∂u/∂z) + ∂/∂z (ν_t ∂w/∂x)
dv/dt = ∂/∂z (ν_t ∂v/∂z) + ∂/∂z (ν_t ∂w/∂y)
dw/dt = ∂/∂z (ν_t ∂w/∂z) - ∂/∂z (ν_t ∂u/∂x) - ∂/∂z (ν_t ∂v/∂y)
```

**Scalar tendencies:**
```cpp
dθ/dt = ∂/∂z (ν_t/Pr_t ∂θ/∂z)
dq/dt = ∂/∂z (ν_t/Pr_t ∂q/∂z)
```

### Boundary Conditions

**Surface fluxes** computed from MOST:
```cpp
τ = -ρ u_*² = ρ ν_t |∂U/∂z|  (momentum flux)
H = -ρ c_p u_* θ_* = ρ c_p ν_t/Pr_t ∂θ/∂z  (heat flux)
E = -ρ u_* q_* = ρ ν_t/Pr_t ∂q/∂z  (moisture flux)
```

**Top boundary**: ν_t = 0 (free slip, zero flux)

### Numerical Stability

**Time step limitations:**
```cpp
Δt_turb ≤ C * (Δz)² / ν_t_max
```

**Implicit treatment** for stiff near-surface gradients:
```cpp
// Crank-Nicolson for vertical diffusion
θ^{n+1} - θ^n = (Δt/2) [ν_t^n ∇²θ^n + ν_t^{n+1} ∇²θ^{n+1}]
```

## Validation and Testing

### Turbulence Metrics
```bash
# Test turbulence schemes
python tests/validate_turbulence.py --scheme tke --resolution 50m

# Compare SGS dissipation rates
python tests/turbulence_diagnostics.py --case neutral_boundary_layer
```

### Expected Behaviors
- **Neutral conditions**: ν_t ∝ z^{4/3} (surface layer scaling)
- **Stable stratification**: Reduced mixing, enhanced gradients
- **Convective conditions**: Enhanced mixing, reduced gradients
- **Energy cascade**: Proper dissipation at SGS scales

## Performance Characteristics

### Computational Cost
- **Smagorinsky**: ~2-5% of total simulation time
- **TKE**: ~5-10% of total simulation time
- **Memory**: Additional 3D arrays for TKE (prognostic schemes)

### Scaling Properties
- **Resolution dependence**: ν_t scales with Δ⁴ (Smagorinsky) or Δ² (TKE)
- **Stability**: TKE schemes more stable for coarse grids
- **Parallel efficiency**: Good scalability for distributed grids

## Scientific References

### Smagorinsky (1963)
- Original eddy-viscosity formulation
- Foundation for atmospheric SGS models
- Lilly (1962) stability corrections

### Deardorff (1980)
- Prognostic TKE for boundary layers
- Extension to convective conditions
- Length scale formulations

### Moeng (1984)
- LES validation of SGS schemes
- Turbulent channel flow benchmarks
- Buoyancy effects on turbulence

### Nakanishi (2001)
- Mellor-Yamada level 2.5 adaptation
- Atmospheric boundary layer optimization
- Stability function improvements

This module provides physically-based sub-grid scale mixing for convective storm simulations.
