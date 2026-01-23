# Numerical Methods Module

This module provides the core numerical algorithms for spatial discretization and time integration in the atmospheric simulation framework.

## Overview

The numerics module implements:
- **Spatial discretization**: Finite-difference and interpolation operators
- **Advection schemes**: Conservative transport of scalar fields
- **Diffusion operators**: Explicit and implicit Laplacian solvers
- **Time integration**: Runge-Kutta methods with stability monitoring
- **Grid operations**: Interpolation and boundary condition handling

## Architecture

```
src/numerics/
├── numerics.cpp              # Module coordinator
├── base/                     # Common utilities (empty)
├── advection/                # Scalar transport schemes
│   ├── factory.cpp/.hpp      # Advection scheme factory
│   └── schemes/
│       ├── tvd/              # Total Variation Diminishing
│       └── weno5/            # Weighted Essentially Non-Oscillatory
├── diffusion/                # Laplacian operators
│   ├── factory.cpp/.hpp      # Diffusion scheme factory
│   └── schemes/
│       ├── explicit/         # Forward Euler diffusion
│       └── implicit/         # Crank-Nicolson implicit
└── time_stepping/            # Time integration methods
    ├── factory.cpp/.hpp      # Time stepping factory
    └── schemes/
        ├── rk3/              # 3rd-order Runge-Kutta
        └── rk4/              # 4th-order Runge-Kutta
```

## Spatial Discretization

### Grid Geometry

**Cylindrical coordinates** (r, θ, z):
```cpp
// Coordinate arrays
r[i] = i * dr + r_min;              // Radial coordinate
theta[j] = j * dtheta;              // Azimuthal coordinate
z[k] = k * dz + z_min;              // Vertical coordinate

// Grid spacings
dr = (r_max - r_min) / (NR - 1);    // Radial spacing
dtheta = 2 * M_PI / NTH;            // Azimuthal spacing
dz = (z_max - z_min) / (NZ - 1);    // Vertical spacing
```

**Staggering**: Lorenz grid (thermodynamics at cell centers, momentum at faces)

### Differential Operators

**Gradient operators:**
```cpp
// Central differences
∂ψ/∂r |_{i,j,k} = (ψ_{i+1,j,k} - ψ_{i-1,j,k}) / (2 dr)
∂ψ/∂θ |_{i,j,k} = (ψ_{i,j+1,k} - ψ_{i,j-1,k}) / (2 r dtheta)
∂ψ/∂z |_{i,j,k} = (ψ_{i,j,k+1} - ψ_{i,j,k-1}) / (2 dz)
```

**Laplacian operator:**
```cpp
∇²ψ |_{i,j,k} = ∂²ψ/∂r² + (1/r)∂ψ/∂r + (1/r²)∂²ψ/∂θ² + ∂²ψ/∂z²
```

## Advection Schemes

### Total Variation Diminishing (TVD)

**Monotonic upstream-centered scheme for conservation laws** (van Leer, 1977):

**1D advection equation:**
```cpp
∂ψ/∂t + u ∂ψ/∂x = 0
```

**TVD discretization:**
```cpp
ψ_i^{n+1} = ψ_i^n - (Δt/Δx) [F_{i+1/2} - F_{i-1/2}]
```

Where F_{i+1/2} is the numerical flux using limiter functions.

**Limiter functions:**
- **minmod**: ϕ(r) = max(0, min(1, r))
- **van Leer**: ϕ(r) = (r + |r|) / (1 + |r|)
- **superbee**: ϕ(r) = max(0, min(2r, 1), min(r, 2))

### Weighted Essentially Non-Oscillatory (WENO5)

**5th-order accurate scheme** for smooth regions with shock capturing:

**WENO reconstruction** (Jiang & Shu, 1996):
```cpp
// Candidate stencils
q^{(0)} = (2/6)ψ_{i-2} - (7/6)ψ_{i-1} + (11/6)ψ_i
q^{(1)} = (-1/6)ψ_{i-1} + (5/6)ψ_i + (2/6)ψ_{i+1}
q^{(2)} = (2/6)ψ_i + (5/6)ψ_{i+1} - (1/6)ψ_{i+2}

// Nonlinear weights
ω_k = α_k / ∑ α_m,    α_k = γ_k / (ε + β_k)²

// Final reconstruction
ψ_{i+1/2} = ∑ ω_k q^{(k)}
```

**Smoothness indicators** β_k detect discontinuities and reduce weights accordingly.

## Diffusion Schemes

### Explicit Laplacian

**Forward Euler time integration:**
```cpp
ψ^{n+1} = ψ^n + Δt κ ∇²ψ^n
```

**Stability condition:** Δt ≤ (Δx²)/(2d κ) where d is dimensionality

### Implicit Crank-Nicolson

**Second-order accurate implicit scheme:**
```cpp
ψ^{n+1} - ψ^n = (Δt/2) κ [∇²ψ^{n+1} + ∇²ψ^n]
```

**Resulting system:** A ψ^{n+1} = ψ^n where A = I - (Δt/2)κ ∇²

**Solution method:** Successive over-relaxation (SOR) or conjugate gradient

## Time Integration

### Runge-Kutta Methods

**3rd-order TVD RK3** (Shu & Osher, 1988):
```cpp
u^{(1)} = u^n + Δt L(u^n)
u^{(2)} = (3/4)u^n + (1/4)u^{(1)} + (1/4)Δt L(u^{(1)})
u^{n+1} = (1/3)u^n + (2/3)u^{(2)} + (2/3)Δt L(u^{(2)})
```

**4th-order classical RK4:**
```cpp
k₁ = Δt L(u^n)
k₂ = Δt L(u^n + k₁/2)
k₃ = Δt L(u^n + k₂/2)
k₄ = Δt L(u^n + k₃)
u^{n+1} = u^n + (k₁ + 2k₂ + 2k₃ + k₄)/6
```

### CFL Monitoring

**Courant-Friedrichs-Lewy condition:**
```cpp
CFL = Δt * max(|u|/Δx + |v|/Δy + |w|/Δz) ≤ 1.0
```

**Adaptive time stepping:**
```cpp
Δt_new = min(Δt_max, CFL_target * Δx / max_speed)
```

## Configuration

### Advection Settings
```yaml
numerics:
  advection:
    scheme: "weno5"       # tvd, weno5
    limiter: "van_leer"   # minmod, van_leer, superbee (TVD only)
    order: 5              # 3 or 5 (WENO)
```

### Diffusion Settings
```yaml
numerics:
  diffusion:
    scheme: "explicit"    # explicit, implicit
    kappa_max: 100.0      # Maximum diffusion coefficient
    tolerance: 1e-6       # Implicit solver tolerance
    max_iterations: 100   # Maximum solver iterations
```

### Time Stepping Settings
```yaml
numerics:
  time_stepping:
    scheme: "rk3"         # rk3, rk4
    cfl_target: 0.5       # Target CFL number
    dt_min: 0.01          # Minimum time step (s)
    dt_max: 10.0          # Maximum time step (s)
```

## Implementation Details

### Boundary Conditions

**Radial boundaries:**
```cpp
// Reflective (walls)
ψ[0] = ψ[1]      // Zero-gradient
u[0] = -u[1]     // Reflective

// Periodic (azimuthal)
ψ[NTH] = ψ[0]    // Continuity
```

**Vertical boundaries:**
```cpp
// Bottom: free-slip
w[0] = 0
∂u/∂z|₀ = 0

// Top: Rayleigh damping
τ(z) = τ₀ exp(-((z-z_damp)/H_damp)²)
```

### Parallel Considerations

**Domain decomposition:**
- Radial direction: Natural for cylindrical geometry
- Azimuthal direction: Periodic boundaries
- Vertical direction: Surface/top boundary handling

**Ghost cell exchange:**
- MPI communication for multi-domain simulations
- Halo updates for finite-difference stencils

### Performance Optimization

**Cache-efficient data layouts:**
```cpp
// Array ordering: field[r][theta][z]
std::vector<std::vector<std::vector<float>>> psi(NR,
    std::vector<std::vector<float>>(NTH,
        std::vector<float>(NZ)));
```

**SIMD vectorization:**
- Inner loops optimized for AVX/AVX-512
- Structure-of-arrays for vector operations

## Validation and Testing

### Numerical Accuracy Tests
```bash
# Test advection schemes
python tests/test_advection.py --scheme weno5 --case solid_body_rotation

# Validate diffusion operators
python tests/test_diffusion.py --scheme implicit --grid 64x64x64
```

### Convergence Studies
- **Order verification**: Error vs. resolution scaling
- **Conservation properties**: Mass, momentum, energy conservation
- **Stability analysis**: Eigenvalue analysis of linearized operators

### Benchmark Results
- **TVD scheme**: 2nd-order accurate, monotonic
- **WENO5 scheme**: 5th-order accurate in smooth regions
- **RK3**: 3rd-order time accuracy, TVD stable
- **Implicit diffusion**: Unconditionally stable, 2nd-order accurate

## Scientific References

### Advection Schemes
- Godunov (1959) - Finite-difference methods for shock waves
- van Leer (1979) - Towards the ultimate conservative difference scheme
- Jiang & Shu (1996) - Efficient implementation of WENO schemes

### Time Integration
- Shu & Osher (1988) - Efficient implementation of essentially non-oscillatory schemes
- Williamson (1980) - Low-storage Runge-Kutta schemes

### Numerical Methods
- Durran (2010) - Numerical Methods for Fluid Dynamics
- Wicker & Skamarock (2002) - Time-splitting methods for compressible models

This module provides the numerical foundation for accurate and stable atmospheric simulations.
