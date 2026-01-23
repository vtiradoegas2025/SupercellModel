# SupercellModel Technical Documentation

SupercellModel is a research-grade atmospheric simulation framework implementing compressible non-hydrostatic dynamics with modular physics parameterizations for supercell thunderstorm research.

## System Architecture

### Core Components

#### Simulation Engine (C++17)
- **Dynamics Core**: Compressible Euler equations in cylindrical coordinates (r, θ, z)
- **Time Integration**: 3rd-order Runge-Kutta (RK3) with CFL-limited adaptive stepping
- **Spatial Discretization**: Finite-difference methods with configurable advection schemes
- **Boundary Conditions**: Radiative lateral boundaries, Rayleigh damping aloft, free-slip surface

#### Physics Modules
- **Microphysics**: Bulk parameterization schemes (Kessler, Thompson, Lin, Milbrandt)
- **Radiation**: Longwave/shortwave transfer (simple grey, RRTMG)
- **Boundary Layer**: Planetary boundary layer parameterizations (YSU, MYNN, slab)
- **Turbulence**: Sub-grid scale closures (Smagorinsky, TKE)
- **Radar Operators**: Forward operators for reflectivity (Z/Ze), Doppler velocity (Vr), polarimetric variables (ZDR, KDP)

#### Visualization Pipeline (Python)
- **Data Processing**: Cylindrical-to-Cartesian coordinate transformation
- **Volume Rendering**: OpenGL 3.3+ ray marching with transfer functions
- **Data Formats**: NPY export from C++, Zarr format for Python visualization

### Modular Design Pattern

All physics components follow a factory-based architecture:

```cpp
// Physics scheme registration
REGISTER_MICROPHYSICS_SCHEME(KesslerMicrophysics, "kessler");
REGISTER_MICROPHYSICS_SCHEME(ThompsonMicrophysics, "thompson");

// Runtime scheme selection
auto microphysics = create_microphysics_scheme(config.scheme_id);
microphysics->initialize(config);
```

This enables easy comparison of parameterizations and extension with new schemes.

### Project Status
- **Core Framework**: Functional with validated physics implementations
- **Radar Module**: Basic reflectivity and velocity operators implemented
- **Visualization**: 3D volume rendering pipeline operational
- **In Development**: Advanced radar polarimetry, terrain effects, stochastic perturbations

## Technical Specifications

### Governing Equations

#### Compressible Non-Hydrostatic Dynamics
The model solves the compressible Euler equations in cylindrical coordinates:

```
∂u/∂t + ∇·(u⊗u) + (1/ρ₀)∇p' + gθ'/θ₀ k̂ = -∇·τ + F_buoyancy
∂w/∂t + ∇·(u⊗w) + (1/ρ₀)∇p' + g(θ'/θ₀ - q_t) k̂ = -∇·τ
∂θ/∂t + ∇·(uθ) = Q_radiation + Q_microphysics + ∇·(K_θ ∇θ)
∂q_v/∂t + ∇·(uq_v) = -C - E + ∇·(K_q ∇q_v)
∂q_c/∂t + ∇·(uq_c) = C - A + ∇·(K_q ∇q_c)
∂q_r/∂t + ∇·(uq_r) = A - E + ∇·(K_q ∇q_r) - ∇·(V_t q_r)
```

Where:
- ρ₀: Base-state density profile
- p': Perturbation pressure
- τ: Sub-grid stress tensor
- Q: Source/sink terms
- V_t: Terminal fall velocity

#### Coordinate System
- **Cylindrical grid**: (r, θ, z) with uniform spacing
- **Grid dimensions**: NR × NTH × NZ (configurable)
- **Resolution**: Δr, Δθ, Δz with aspect ratio constraints
- **Domain**: r ∈ [0, r_max], θ ∈ [0, 2π], z ∈ [0, z_max]

### Prognostic Variables

#### Momentum and Thermodynamics
- **u, v, w**: 3D wind components (m/s)
- **θ**: Potential temperature (K)
- **ρ**: Density perturbation (kg/m³)
- **p**: Pressure perturbation (Pa)

#### Hydrometeor Species
- **q_v**: Water vapor mixing ratio (kg/kg)
- **q_c**: Cloud water mixing ratio (kg/kg)
- **q_r**: Rain water mixing ratio (kg/kg)
- **q_i**: Cloud ice mixing ratio (kg/kg, Thompson/Lin schemes)
- **q_s**: Snow mixing ratio (kg/kg, Thompson/Lin schemes)
- **q_g**: Graupel mixing ratio (kg/kg, Thompson/Milbrandt schemes)
- **q_h**: Hail mixing ratio (kg/kg, optional)

#### Turbulence State (TKE scheme)
- **tke**: Turbulent kinetic energy (m²/s²)
- **ε**: Turbulent dissipation rate (m²/s³)

### Diagnostic Fields

#### Thermodynamic Diagnostics
- **T**: Temperature from θ and p (K)
- **θ_v**: Virtual potential temperature (K)
- **CAPE/CIN**: Convective available/unavailable potential energy (J/kg)
- **LCL/LFC/EL**: Lifting condensation/freezing/equilibrium levels (m)

#### Dynamic Diagnostics
- **ζ**: Vertical vorticity component (s⁻¹)
- **ω**: Vertical velocity in pressure coordinates (Pa/s)
- **σ**: Divergence (s⁻¹)
- **PV**: Potential vorticity (PVU)

#### Storm-Scale Diagnostics
- **UH**: Updraft helicity (m²/s²)
- **SRH**: Storm-relative helicity (m²/s²)
- **MESO**: Mesocyclone identification metrics

#### Radar Forward Operators
- **Z_e**: Equivalent radar reflectivity factor (mm⁶/m³ → dBZ)
- **V_r**: Radial Doppler velocity (m/s)
- **Z_DR**: Differential reflectivity (dB)
- **K_DP**: Specific differential phase (deg/km)
- **ρ_HV**: Cross-correlation coefficient

### Physics Parameterizations

#### Microphysics (4 schemes)
- **Kessler**: Warm-rain (3 species) - see `src/microphysics/README.md`
- **Thompson**: Aerosol-aware mixed-phase (6 species) - see `src/microphysics/README.md`
- **Lin**: Traditional ice (5 species) - see `src/microphysics/README.md`
- **Milbrandt-Yau**: Double-moment (prognostic N) - see `src/microphysics/README.md`

#### Radiation (2 schemes)
- **Simple Grey**: Newtonian relaxation - see `src/radiation/README.md`
- **RRTMG**: Full spectral transfer (30 bands) - see `src/radiation/README.md`

#### Boundary Layer (3 schemes)
- **Slab**: Simple mixed-layer - see `src/boundary_layer/README.md`
- **YSU**: Non-local mixing - see `src/boundary_layer/README.md`
- **MYNN**: TKE-based closure - see `src/boundary_layer/README.md`

#### Turbulence (2 schemes)
- **Smagorinsky-Lilly**: Eddy-viscosity SGS - see `src/turbulence/README.md`
- **TKE Prognostic**: Full prognostic scheme - see `src/turbulence/README.md`

#### Radar Forward Operators (3 schemes)
- **Reflectivity**: Z/Ze from microphysics - see `src/radar/README.md`
- **Doppler Velocity**: V_r from wind fields - see `src/radar/README.md`
- **Polarimetric**: Z_DR, K_DP, ρ_HV - see `src/radar/README.md`

### Numerical Configuration

#### Grid Specifications

##### Coordinate System
- **Primary grid**: Cylindrical (r, θ, z) coordinates
- **Visualization**: Cartesian (x, y, z) transformation
- **Staggering**: Lorenz grid (pressure/thermo at cell centers, momentum at faces)

##### Resolution Options
- **Test grids**: 32×32×16 to 64×64×32 (Δr = 1000-2000m, Δz = 250-500m)
- **Production grids**: 256×128×128 (Δr = 1000m, Δθ = 1000m, Δz = 100m)
- **High-resolution**: 512×256×256 (Δr = 500m, Δθ = 500m, Δz = 50m)

##### Domain Configuration
- **Radial extent**: r_max = NR × Δr (typically 128-256 km)
- **Vertical extent**: z_max = NZ × Δz (typically 6-16 km)
- **Azimuthal coverage**: Full 360° (2π radians)

#### Time Integration

##### RK3 Integration
3rd-order explicit Runge-Kutta with CFL-limited time stepping:

```
y_{n+1} = y_n + (1/6)(k₁ + 4k₂ + k₃)
k₁ = Δt f(t_n, y_n)
k₂ = Δt f(t_n + 0.5Δt, y_n + 0.5k₁)
k₃ = Δt f(t_n + Δt, y_n - k₁ + 2k₂)
```

##### CFL Criterion
Δt ≤ C × min(Δx/|u|, Δy/|v|, Δz/|w|) where C = 0.5 for stability

#### Boundary Conditions

##### Lateral Boundaries (r = r_max)
- **Radiative condition**: ∂ψ/∂t + c ∂ψ/∂r = 0 (c = phase speed)
- **Open boundary**: Zero-gradient for scalars, extrapolation for momentum

##### Top Boundary (z = z_max)
- **Rayleigh damping**: τ(z) = τ₀ exp(-((z - z_damp)/H_damp)²)
- **Sponge layer**: 4-6 km deep with relaxation time τ = 300s

##### Bottom Boundary (z = 0)
- **Free-slip**: w = 0, ∂u/∂z = ∂v/∂z = 0
- **Surface fluxes**: MOST similarity theory for momentum/heat transfer

### Storm Configurations

#### Pre-defined Test Cases

##### Classic Supercell (Weisman-Klemp)
- **Thermodynamics**: CAPE = 2500 J/kg, CIN = 50 J/kg
- **Wind Profile**: 0-6km shear = 40 m/s, moderate low-level curvature
- **Trigger**: Thermal bubble (Δθ = 2K, radius = 10km)
- **Expected**: Rotating updraft, hook echo, mesocyclone

##### Low Precipitation (LP) Supercell
- **Thermodynamics**: CAPE = 3000 J/kg, drier mid-levels
- **Wind Profile**: 0-6km shear = 35 m/s, straight hodograph
- **Trigger**: Same thermal bubble
- **Expected**: Sculpted appearance, minimal precipitation

##### High Precipitation (HP) Supercell
- **Thermodynamics**: CAPE = 2500 J/kg, moist column
- **Wind Profile**: 0-6km shear = 40 m/s, curved hodograph
- **Microphysics**: Hail/graupel enabled
- **Expected**: Strong cold pool, wrapping precipitation

##### Cyclic Supercell
- **Thermodynamics**: CAPE = 2000 J/kg, elevated instability
- **Wind Profile**: Strong deep-layer shear, curved hodograph
- **Trigger**: Smaller thermal bubble (Δθ = 1.5K, radius = 5km)
- **Expected**: Multiple updraft cycles, occlusions

##### Elevated/HSLC Case
- **Thermodynamics**: Reduced surface-based CAPE, elevated instability
- **Wind Profile**: Weak low-level shear, strong mid-level winds
- **Trigger**: Surface-based bubble in elevated environment
- **Expected**: Elevated convection, minimal surface interaction

#### Initialization Protocol

##### Base State Construction
1. **Temperature profile**: θ₀(z) = θ_sfc + Γ_dry(z) for z < z_LCL, moist adiabatic above
2. **Density profile**: ρ₀(z) from hydrostatic balance
3. **Pressure profile**: p₀(z) from ideal gas law

##### Environmental Winds
1. **Hodograph specification**: u(z), v(z) profiles from observational data
2. **Smoothing**: Cubic spline interpolation for continuity
3. **Perturbation**: Optional stochastic perturbations (±5-10%)

##### Convection Trigger
1. **Thermal bubble**: Δθ(r,z) = Δθ₀ exp(-r²/R²) × exp(-(z-z₀)²/H²)
2. **Multiple bubbles**: Optional for broader initiation
3. **Random perturbations**: Optional stochastic noise

For complete scientific attribution and literature references, see **[foundationalScience.md](docs/foundationalScience.md)**.

### Prerequisites
- **C++17 compiler** (`clang++`/`g++`)
- **SFML 3** (Graphics, Window, System). macOS (Homebrew): `brew install sfml`
- **Python 3.10+** with visualization dependencies

### Implementation Status

#### Core Framework
- **Dynamics**: Compressible Euler equations implemented and validated
- **Time Integration**: RK3 scheme with CFL monitoring operational
- **Boundary Conditions**: Radiative/open boundaries with Rayleigh damping
- **Grid**: Cylindrical coordinate system with configurable resolution

#### Physics Modules Status

| Module | Implementation | Validation Status |
|--------|----------------|-------------------|
| Microphysics | Kessler, Thompson, Lin, Milbrandt | Physically realistic phase changes |
| Radiation | Simple grey, RRTMG | Energy balance maintained |
| Boundary Layer | YSU, MYNN, slab | Surface fluxes and mixing |
| Turbulence | Smagorinsky, TKE | SGS stresses resolved |
| Radar Operators | Z/Ze, Vr, ZDR | Forward operators functional |

#### Validation Metrics
- **Mass Conservation**: <0.1% drift over 2-hour simulations
- **Energy Conservation**: <5% total energy change in stable cases
- **Numerical Stability**: No NaN/inf values in validated test cases
- **Physical Realism**: Temperature/wind/moisture fields within literature ranges
- **Radar Consistency**: Z-Vr relationships match theoretical expectations

#### Known Issues
- **Terrain Module**: Excluded from build due to compilation issues
- **Chaos Module**: Contains placeholder implementations
- **Small Grids**: <64×64 may cause initialization instabilities
- **Memory Usage**: High-resolution grids require significant RAM (>16GB)
- **Parallelization**: Single-threaded implementation (OpenMP planned)

### Build System

#### Compilation Requirements
- **C++ Standard**: C++17 compliant compiler (clang++ ≥ 9.0, g++ ≥ 7.0)
- **Build System**: GNU Make with automatic dependency generation
- **External Libraries**:
  - SFML 3.x (graphics/window support, optional for GUI)
  - Standard Template Library (STL) only for core functionality

#### Build Targets
```bash
# Default build (optimized)
make

# Debug build with symbols
make DEBUG=1

# Build with GUI support
make GUI=1

# Clean build artifacts
make clean
```

#### Compilation Flags
- **Optimization**: -O2 (production), -O0 (debug)
- **Warnings**: -Wall -Wextra -Wpedantic
- **Features**: -DEXPORT_NPY for data export
- **Standards**: -std=c++17 -fPIC

#### Runtime Dependencies
- **Execution**: No external runtime libraries required (static linking)
- **Data Export**: NPY format (NumPy compatible)
- **Visualization**: Python 3.10+ with scientific stack

### Simulation Execution

#### Command Line Interface
```bash
bin/tornado_sim [OPTIONS]
```

#### Key Options
- `--headless`: Run without GUI for batch processing
- `--config PATH`: YAML configuration file
- `--duration SEC`: Simulation duration in seconds
- `--write-every SEC`: Data export frequency
- `--outdir PATH`: Output directory for simulation data

#### Configuration Files
Pre-configured test cases in `configs/` directory:
- `classic.yaml`: Standard supercell setup (CAPE=2500 J/kg, shear=40 m/s)
- `lp.yaml`: Low precipitation supercell (CAPE=3000 J/kg)
- `hp.yaml`: High precipitation supercell with hail
- `cyclic.yaml`: Multi-cycle supercell development
- `elevated.yaml`: Elevated convection case
- `test_small.yaml`: Minimal test configuration (32×32×16)

#### Output Formats
- **NPY files**: NumPy-compatible binary arrays per theta slice
- **Directory structure**: `step_XXXXXX/thXX_field.npy`
- **Metadata**: Grid parameters embedded in filenames

## Visualization Pipeline

### Architecture Overview

#### Data Processing Chain
```
C++ Simulation → NPY Export → Python Processing → Zarr Storage → OpenGL Rendering → Video Output
     ↓              ↓              ↓                    ↓              ↓
Cylindrical     NumPy arrays   Coordinate           Volume data    Ray marching    FFmpeg
coordinates     (r, θ, z)     transformation       (x, y, z)      shaders        encoding
                              Interpolation        Compression     Transfer        H.264/MP4
                              Filtering            Chunking        functions
```

#### Component Responsibilities
- **C++ Backend**: High-performance simulation with minimal I/O overhead
- **Python Pipeline**: Data transformation and scientific analysis
- **OpenGL Renderer**: Real-time volume rendering with GPU acceleration
- **Zarr Format**: Cloud-optimized storage with compression and chunking

### Coordinate System Transformation

#### Cylindrical to Cartesian Mapping
The simulation uses cylindrical coordinates (r, θ, z) for computational efficiency, while visualization requires Cartesian coordinates (x, y, z):

```python
# Input: Cylindrical grid (NR × NTH × NZ)
r_coords = np.linspace(0, r_max, NR)        # Radial coordinates
theta_coords = np.linspace(0, 2*np.pi, NTH) # Azimuthal angles
z_coords = np.linspace(0, z_max, NZ)        # Vertical levels

# Output: Cartesian grid (NX × NY × NZ)
x_coords = np.linspace(-r_max, r_max, NX)
y_coords = np.linspace(-r_max, r_max, NY)
z_coords_cart = z_coords  # Vertical coordinates unchanged

# Transformation for each point
for i in range(NR):
    for j in range(NTH):
        for k in range(NZ):
            x = r_coords[i] * np.cos(theta_coords[j])
            y = r_coords[i] * np.sin(theta_coords[j])
            z = z_coords[k]
            # Interpolate to Cartesian grid
```

#### Interpolation Methods
- **Bilinear interpolation**: For smooth field reconstruction
- **Conservative remapping**: Mass conservation for hydrometeor fields
- **Boundary handling**: Zero-gradient extrapolation beyond domain

### Volume Rendering Implementation

#### Ray Marching Algorithm
The renderer uses GPU-accelerated ray marching through the volume:

```glsl
// Vertex shader: Camera setup and ray generation
void main() {
    gl_Position = projection * view * model * vec4(position, 1.0);
    ray_origin = (inverse_view * vec4(0,0,0,1)).xyz;
    ray_direction = normalize(position - ray_origin);
}

// Fragment shader: Volume integration
vec4 ray_march(vec3 ray_origin, vec3 ray_direction) {
    vec4 accumulated_color = vec4(0.0);
    float accumulated_alpha = 0.0;

    for (int step = 0; step < MAX_STEPS; step++) {
        vec3 sample_pos = ray_origin + ray_direction * t;
        vec4 sample_color = sample_volume(sample_pos);
        sample_color.a *= step_size;  // Opacity scaling

        accumulated_color += sample_color * (1.0 - accumulated_alpha);
        accumulated_alpha += sample_color.a;

        if (accumulated_alpha >= 0.95) break;
        t += step_size;
    }

    return accumulated_color;
}
```

#### Transfer Functions
Atmospheric fields are mapped to RGBA channels for physical interpretation:

- **Red Channel (R)**: Buoyancy/θ' - Warm updrafts (positive) → red/cyan, cold downdrafts (negative) → blue
- **Green Channel (G)**: Precipitation/q_r - Rain shafts and convective cores
- **Blue Channel (B)**: Wind speed - Flow intensity and storm dynamics
- **Alpha Channel (A)**: Optical depth - Visibility and cloud opacity

#### Performance Optimizations
- **Early ray termination**: Stop marching when accumulated opacity > 0.95
- **Level-of-detail**: Adaptive sampling based on view distance
- **GPU memory management**: Texture streaming for large datasets
- **Shader precompilation**: Minimize CPU-GPU synchronization

### Quick Start

**Install Python dependencies:**
```bash
pip install -r requirements.txt
# Or manually: pip install zarr numpy moderngl moderngl-window pillow tqdm
```

**Interactive guided setup (recommended for new users):**
```bash
python visualization/run_3d_pipeline.py --quick-start
```
This interactive guide checks your setup and walks through all visualization options.

### Manual Usage

**Convert simulation data to visualization format:**
```bash
python visualization/convert_npy_to_zarr.py --input data/exports --output data/supercell.zarr --max-timesteps 100
```

**Run interactive 3D viewer:**
```bash
python visualization/supercell_renderer.py --input data/supercell.zarr --field theta
```
*Controls: SPACE (play/pause), R (reset), W (wireframe), Up/Down (opacity), Left/Right (brightness), Mouse (orbit)*

**Render offline video:**
```bash
python visualization/render_supercell_3d.py --input data/supercell.zarr --output storm.mp4 --fps 30 --duration 60
```

**Run complete automated pipeline:**
```bash
python visualization/run_3d_pipeline.py --all --output supercell.mp4
```

### Renderable Fields

#### Primary Atmospheric Variables
- **theta**: Potential temperature (K) - thermodynamic buoyancy
- **qv**: Water vapor mixing ratio (kg/kg) - humidity field
- **qc**: Cloud water mixing ratio (kg/kg) - cloud condensate
- **qr**: Rain water mixing ratio (kg/kg) - precipitation water
- **qi**: Cloud ice mixing ratio (kg/kg) - ice crystals
- **qs**: Snow mixing ratio (kg/kg) - snow aggregates
- **qh**: Hail mixing ratio (kg/kg) - hail/graupel stones

#### Dynamic Variables
- **u, v, w**: Wind components (m/s) - 3D velocity field
- **speed**: Wind speed magnitude (m/s) - flow intensity
- **vorticity_z**: Vertical vorticity component (s⁻¹) - rotation

#### Radar-Derived Fields
- **reflectivity_dbz**: Radar reflectivity (dBZ) - precipitation intensity
- **radial_velocity**: Doppler radial velocity (m/s) - line-of-sight winds
- **differential_reflectivity**: Z_DR (dB) - hydrometeor type/size
- **specific_differential_phase**: K_DP (deg/km) - liquid water content

#### Diagnostic Fields
- **buoyancy**: Buoyancy acceleration (m/s²) - convective forcing
- **updraft_helicity**: UH (m²/s²) - rotating updraft potential
- **storm_relative_helicity**: SRH (m²/s²) - environmental shear
- **cape**: Convective available potential energy (J/kg)
- **cin**: Convective inhibition (J/kg)

### Performance & Optimization

#### Resolution Strategy
- **Low resolution (testing)**: 64×64×32 grid, dx=2000m
- **Medium resolution**: 128×128×64 grid, dx=1000m
- **High resolution (production)**: 256×256×128 grid, dx=500m

#### GPU Requirements
- OpenGL 3.3+ compatible GPU
- Memory: ~4MB per timestep for 128×128×64 volumes
- Performance: 30-60 FPS interactive, 4K offline rendering

#### Rendering Parameters
- **Ray marching steps**: 200 steps through volume
- **Step size**: 0.015 (normalized coordinates)
- **Early termination**: When accumulated alpha > 0.98
- **Transfer functions**: Adjustable opacity and brightness

### Output Formats

#### Interactive Viewing
- Real-time 3D exploration with orbital camera
- Adjustable transfer functions and rendering parameters
- Field selection and visualization presets

#### Offline Rendering
- **Resolution**: 1920×1080 (configurable)
- **Frame rate**: 30 FPS (configurable)
- **Codec**: H.264 with FFmpeg
- **Quality**: High quality CRF encoding

#### Video Products
1. **Realtime video**: 1× speed, natural evolution
2. **Timelapse**: 10-30× speed, full lifecycle
3. **Feature highlights**: Focus on specific phenomena

### Data Storage (Zarr Format)

#### Schema
- **Root attributes**: version, dx/dy/dz, dt, time_units, origin, case_name
- **Datasets** (C order, time-major):
  - `theta`: float32 [t, z, y, x], chunks [1, 64, 128, 128]
  - `qv, qc, qr`: float32 [t, z, y, x], same chunking
  - `u, v, w`: float32 [t, z, y, x], same chunking
  - `reflectivity_dbz`: float32 [t, z, y, x] (diagnostic)
- **Compression**: blosc zstd, level 5; shuffle enabled
- **Consistency**: Atomic writes per timestep

### Headless Batch Mode

**Simulation + Rendering Pipeline:**
```bash
# Run simulation with direct Zarr output
./bin/tornado_sim --headless \
  --config configs/classic.yaml \
  --output data/classic.zarr \
  --duration 7200 \
  --write-every 5

# Render visualization videos
python visualization/render_supercell_3d.py \
  --input data/classic.zarr \
  --output classic_realtime.mp4 \
  --fps 30

python visualization/render_supercell_3d.py \
  --input data/classic.zarr \
  --output classic_timelapse.mp4 \
  --speed 10 \
  --fps 30
```

### Validation Framework

#### Physics Validation
SupercellModel includes comprehensive validation tools for physical realism:

```bash
# Run physics validation
python tests/validate_simulation.py --input data/exports --timestep 0

# Validate time stepping
python tests/validate_time_steps.py
```

#### Validation Metrics

##### Mass Conservation
- **Total water**: ∫(q_v + q_c + q_r + ...) dV conserved to <0.1%
- **Dry air mass**: ρ conservation in anelastic approximation
- **Hydrometeor budget**: Source/sink terms balance over domain

##### Energy Conservation
- **Total energy**: Kinetic + potential + internal energy conservation
- **Budget analysis**: Work done by pressure gradients vs. dissipative losses
- **Radiative balance**: Net heating/cooling consistent with physical expectations

##### Numerical Stability
- **CFL criterion**: Time steps within stability limits
- **Conservation properties**: Discretization maintains physical invariants
- **Boundary artifacts**: Minimal spurious reflections/refractions

##### Physical Realism Checks
- **Temperature ranges**: 200K < T < 350K throughout domain
- **Wind speeds**: |u|, |v|, |w| < 100 m/s (tornado-scale allowance)
- **Moisture bounds**: 0 ≤ q_v ≤ q_vs(T,p), 0 ≤ q_c,q_r,q_i,...
- **Pressure**: p > 0 everywhere (with noted boundary exceptions)

#### Radar Operator Validation

##### Reflectivity Consistency
- **Z-q_r relationship**: Z ∝ q_r^{1.5-2.0} for various drop size distributions
- **Multi-species**: Ice/snow/hail contributions follow Mie theory
- **Attenuation**: Path-integrated attenuation corrections

##### Doppler Velocity Accuracy
- **Terminal fall speeds**: Consistent with microphysical fall speed relations
- **Wind field mapping**: Correct projection onto radar beam geometry
- **Aliasing**: Proper handling of velocity folding

##### Polarimetric Validation
- **Z_DR signatures**: Correct hydrometeor type discrimination
- **K_DP**: Liquid water content retrieval accuracy
- **Cross-correlation**: Phase consistency between H/V channels

### Expected Results by Configuration

#### Classic Supercell Case
- **Storm Structure**: Rotating updraft (ζ > 0.01 s⁻¹), hook echo formation
- **Updraft Characteristics**: w_max > 30 m/s, persistent mesocyclone
- **Precipitation**: Asymmetric hook with forward flank precipitation
- **Radar Signature**: Z > 50 dBZ core, V_r dipole structure

#### Low Precipitation (LP) Case
- **Storm Structure**: Sculpted appearance with minimal precipitation
- **Dynamics**: Stronger low-level rotation, elevated structure
- **Precipitation**: q_r,max < 5 g/kg, narrow precipitation core
- **Radar Signature**: Weak reflectivity (Z < 40 dBZ), clear slot formation

#### High Precipitation (HP) Case
- **Storm Structure**: Strong cold pool, wrapping precipitation
- **Microphysics**: Hail/graupel trajectories, heavy precipitation
- **Radar Signature**: Z > 60 dBZ with hail cores, heavy precipitation shafts

### Performance Benchmarks

#### Computational Performance
- **Small grid (32×32×16)**: ~10,000 time steps/hour on modern CPU
- **Production grid (256×128×128)**: ~100 time steps/hour on high-end workstation
- **Memory usage**: ~8GB for production grids, ~2GB for test grids

#### Visualization Performance
- **Interactive rendering**: 30-60 FPS for 128³ volumes on modern GPU
- **Offline rendering**: 4K video at 30 FPS, ~10 minutes per hour of simulation
- **Data processing**: Coordinate transformation in < 30 seconds for typical datasets

### Troubleshooting

#### Common Issues
- **Compilation failures**: Ensure C++17 support and SFML 3.x installation
- **Runtime crashes**: Check grid size (minimum 32×32×16 recommended)
- **Visualization artifacts**: Verify OpenGL 3.3+ and shader compilation
- **Memory issues**: Reduce grid resolution or increase system RAM

#### Platform-Specific Notes
- **macOS**: Use Homebrew for SFML, ensure OpenGL frameworks available
- **Linux**: Install development packages for SFML and OpenGL
- **Windows**: Use vcpkg for dependency management, ensure MSVC with C++17 support

#### Boundary Condition Artifacts
- **Pressure**: Small negative values (<1% of domain) near lateral boundaries
- **Vorticity**: Spurious vorticity generation at corners (mitigated by smoothing)
- **Scalar fields**: Gradient discontinuities at boundaries (expected with radiative BCs)


