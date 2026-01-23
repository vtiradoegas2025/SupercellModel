# Terrain/Orography Module

This module handles topographic effects on atmospheric flow, including idealized mountain profiles and coordinate transformations for terrain-following grids.

## Overview

Terrain effects are crucial in atmospheric simulations as they can:
- Generate mountain waves and gravity waves
- Modify boundary layer flows
- Influence precipitation patterns
- Create flow blocking and deflection

## Architecture

```
src/terrain/
├── terrain.cpp               # Module coordinator
├── factory.cpp/.hpp          # Terrain scheme factory
├── base/
│   └── topography.cpp/.hpp   # Terrain generation utilities
└── schemes/                  # Terrain implementations
    ├── none/                 # Flat terrain
    ├── bell/                 # Axisymmetric bell mountain
    └── schar/                # Schar-type ridge
```

## Terrain Representations

### Analytical Profiles

#### Bell Mountain
Classic axisymmetric mountain profile:

```cpp
h(r) = h_0 / (1 + (r/R)^2)^1.5
```

Where:
- h_0: Maximum mountain height
- R: Mountain half-width
- r: Radial distance from mountain center

#### Schar Ridge
2D ridge profile for studying gravity waves:

```cpp
h(x) = h_0 exp(-(x/L)^2) cos²(π x / λ)
```

Where:
- L: Mountain half-width
- λ: Wavelength of small-scale features
- x: Along-ridge coordinate

### Coordinate Transformations

#### Terrain-Following Coordinates
Transform from Cartesian (x, y, z) to terrain-following (x, y, σ):

```cpp
σ = (z - h(x,y)) / (H - h(x,y))
```

Where:
- σ ∈ [0, 1]: Terrain-following coordinate
- H: Model top height
- h(x,y): Terrain height field

#### Metric Terms
Coordinate transformation introduces metric terms:

```cpp
∂/∂z = (1/(H - h)) ∂/∂σ
∂/∂x = ∂/∂x - (∂h/∂x)/(H - h) ∂/∂σ
```

## Scheme Implementations

### None (Flat Terrain)

**Configuration:**
```yaml
terrain:
  scheme: "none"
```

**Features:**
- h(x,y) = 0 everywhere
- Cartesian coordinates
- No metric terms
- Simplest configuration

### Bell Mountain

**Configuration:**
```yaml
terrain:
  scheme: "bell"
  height_m: 1000.0         # Maximum height
  half_width_m: 5000.0     # Mountain half-width
  center_x_km: 50.0        # Center position
  center_y_km: 50.0
```

**Applications:**
- Idealized gravity wave studies
- Mountain wave drag parameterization
- Symmetric flow responses

### Schar Ridge

**Configuration:**
```yaml
terrain:
  scheme: "schar"
  height_m: 500.0          # Maximum height
  half_width_m: 5000.0     # Mountain half-width
  wavelength_m: 4000.0     # Small-scale wavelength
  center_x_km: 25.0        # Center position
```

**Features:**
- 2D ridge geometry
- Small-scale perturbations
- Controlled wavelength for resonance studies

## Implementation Details

### Grid Generation

**Terrain field initialization:**
```cpp
for (int i = 0; i < NR; ++i) {
    for (int j = 0; j < NTH; ++j) {
        double x = i * dr;
        double y = j * dtheta * r_mean;
        h[i][j] = compute_terrain_height(x, y);
    }
}
```

**Coordinate transformation:**
```cpp
// Physical height to terrain-following
sigma[k] = (z[k] - h[i][j]) / (z_max - h[i][j]);

// Terrain-following to physical height
z[k] = h[i][j] + sigma[k] * (z_max - h[i][j]);
```

### Metric Terms Calculation

**Vertical metric:**
```cpp
dz_dsigma = z_max - h[i][j];  // Physical layer thickness
dsigma_dz = 1.0 / dz_dsigma;  // Transformation Jacobian
```

**Horizontal metrics:**
```cpp
// Slope terms
dh_dx = compute_gradient(h, DX, i, j);
dh_dy = compute_gradient(h, DY, i, j);

// Full metric coefficients
G11 = 1.0;  // No stretching in x
G12 = 0.0;  // No rotation
G22 = 1.0;  // No stretching in y
G13 = -dh_dx * dsigma_dz;  // Slope term
G23 = -dh_dy * dsigma_dz;  // Slope term
G33 = dsigma_dz * dsigma_dz;  // Vertical stretching
```

### Momentum Equation Modifications

**Terrain-following momentum equations:**
```cpp
// Horizontal momentum
Du/Dt = -∇_σ p'/ρ + ν ∇_σ² u + F_u

// Vertical momentum (terrain-following)
Dw/Dt = -∂p'/∂σ / (ρ (H-h)) + B + ν ∇_σ² w + F_w
```

Where ∇_σ is the terrain-following gradient operator.

### Boundary Conditions

**Lower boundary:**
```cpp
// Terrain surface
w_σ = 0;  // No flow through terrain

// Free-slip or no-slip
∂u/∂σ = 0;  // Free-slip condition
```

**Upper boundary:**
```cpp
// Model top
w_σ = w_top;  // Specified or zero
```

## Configuration

### Basic Terrain Setup
```yaml
terrain:
  scheme: "bell"           # none, bell, schar
  coordinate_system: "cartesian"  # cartesian, terrain_following
```

### Advanced Options
```yaml
terrain:
  # Coordinate system
  terrain_following: true

  # Numerical options
  smoothing: true          # Smooth terrain gradients
  smoothing_iterations: 5

  # Physical options
  blocking_effect: true    # Include flow blocking
  gravity_wave_drag: false # Explicit drag parameterization
```

## Validation and Testing

### Terrain Profile Tests
```bash
# Test terrain generation
python tests/test_terrain.py --scheme bell --height 1000 --width 5000

# Validate coordinate transformation
python tests/test_coordinates.py --terrain bell --resolution 128x128
```

### Flow Response Validation
- **Linear theory**: Mountain wave amplitudes and wavelengths
- **Nonlinear effects**: Flow splitting and wave breaking
- **Boundary layer**: Modified surface fluxes over complex terrain

### Expected Results

#### Bell Mountain Flow
- **Upwind**: Flow deceleration, speedup over mountain
- **Lee waves**: Gravity wave train downwind
- **Rotor formation**: Low-level reversed flow in valleys

#### Schar Ridge Waves
- **Resonant waves**: Controlled wavelength for maximum amplitude
- **Breaking waves**: Turbulent dissipation in upper atmosphere
- **Momentum flux**: Wave drag on large-scale flow

## Performance Considerations

### Computational Cost
- **Flat terrain**: Minimal overhead (baseline)
- **Simple profiles**: O(N_grid) initialization
- **Complex terrain**: Additional metric term calculations

### Memory Requirements
- **Terrain field**: 2D array (NR × NTH)
- **Metric terms**: 3D arrays for transformation coefficients
- **Additional**: Coordinate transformation arrays

## Scientific References

### Terrain-Following Coordinates
- Gal-Chen & Somerville (1975) - Terrain-following coordinates
- Schar et al. (2002) - Terrain-following vertical coordinate
- Klemp (2011) - Smoothed coordinate surfaces

### Mountain Waves and Dynamics
- Smith (1980) - Linear theory of mountain waves
- Durran (1990) - Mountain Meteorology
- Klemp & Lilly (1975) - Severe downslope windstorms

### Flow Over Complex Terrain
- Jackson & Hunt (1975) - Turbulent wind flow over hills
- Belcher & Hunt (1993) - Turbulent flow over hills and waves
- Wood (2000) - Wind flow over complex terrain

This module provides the foundation for simulating atmospheric flow over complex terrain in research applications.
