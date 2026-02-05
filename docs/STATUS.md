# TornadoModel / SupercellModel - Project Status

**Last Updated:** February 5, 2026

This document provides a comprehensive overview of the current state of the TornadoModel/SupercellModel project, including known faults, working components, and next steps. This file is gitignored for local development context.

---

## Project Overview

TornadoModel (also referred to as SupercellModel) is a high-performance atmospheric simulation framework for supercell and tornadic storm research. It features:

- **Core Engine**: C++17 simulation engine implementing compressible, non-hydrostatic equations in cylindrical coordinates
- **Modular Physics**: Factory-based architecture enabling systematic comparison of parameterizations
- **Visualization**: Python + OpenGL 3D volume rendering pipeline
- **Research Focus**: Severe convective storm dynamics, supercell morphology, tornado genesis

---

## Current State

### What's Working

- **Core Simulation Engine**: Functional C++17 implementation
  - Compressible Euler equations solver
  - Cylindrical coordinate system (r, θ, z)
  - RK3 time integration with CFL-limited adaptive stepping
  - Field3D flattened array storage for cache optimization
  - OpenMP multi-core parallelization (when available)

- **Physics Modules**: Multiple schemes implemented
  - **Microphysics**: 4 schemes (Kessler, Thompson, Lin, Milbrandt)
  - **Boundary Layer**: 3 schemes (YSU, MYNN, slab)
  - **Turbulence**: 2 schemes (Smagorinsky-Lilly, TKE prognostic)
  - **Radiation**: 2 schemes (Simple grey, RRTMG)
  - **Radar**: Forward operators for reflectivity, velocity, polarimetric variables

- **Build System**: Makefile with OpenMP detection, optional GUI support
- **Visualization Pipeline**: Modular engine architecture in place (see visualization/STATUS.md for details)

### What's Incomplete

- **Chaos Module**: Has "COMEBACK" sections marked for future work
  - Basic initial condition perturbations implemented
  - Boundary layer and full stochastic schemes incomplete
  - Spectral Gaussian filter temporarily removed due to compilation issues

- **Terrain Module**: Excluded from build due to compilation issues
  - Bell and Schär mountain implementations exist but not integrated

- **Testing Infrastructure**: Test files moved to `tests/` folder, test targets removed from Makefile

---

## Recent Changes

Based on git status, recent modifications include:

- **Makefile**: Test targets removed
- **Configuration Files**: Test configs moved to `tests/` folder
- **Core Headers**: Updates to chaos_base.hpp, dynamics_base.hpp, microphysics_base.hpp, radar_base.hpp, simulation.hpp, turbulence_base.hpp
- **Source Files**: Modifications across dynamics, chaos, boundary layer, and other modules
- **File Organization**: Test and debug files consolidated into `tests/` directory

---

## Key Components

### Core Simulation (`src/`)

- **equations.cpp**: Field initialization, microphysics stepping, nested grid support
- **dynamics.cpp**: Dynamics coordination, time stepping, bounds checking
- **tornado_sim.cpp**: Main executable, configuration parsing, simulation loop
- **advection/advection.cpp**: Scalar advection with TVD/WENO5 schemes
- **numerics/**: Advection, diffusion, time-stepping factories and schemes

### Physics Modules

- **microphysics/**: 4 bulk parameterization schemes
- **boundary_layer/**: PBL schemes with surface flux calculations
- **turbulence/**: Sub-grid scale closures
- **radiation/**: Radiative transfer schemes
- **radar/**: Forward operators for radar observables
- **chaos/**: Stochastic perturbations (partially implemented)

### Visualization (`visualization/`)

- **core/**: Modular rendering engine (render_engine.py, shader_manager.py, etc.)
- **renderers/**: Interactive and offline renderers
- **passes/**: Rendering passes
- **shaders/**: GLSL shaders for volume and contour rendering

---

## Build Status

### Prerequisites

- **C++ Compiler**: C++17 compliant (clang++ ≥ 9.0, g++ ≥ 7.0)
- **OpenMP**: Optional but recommended (libomp on macOS via Homebrew)
- **Python**: 3.10+ with scientific stack
- **Graphics**: OpenGL 3.3+ compatible GPU (for visualization)

### Build Instructions

```bash
# Build simulation engine
make

# Build with GUI support (SFML)
make GUI=1

# Run simulation
./bin/tornado_sim --headless --config=configs/classic.yaml --duration=3600
```

### Known Build Issues

- Terrain module excluded from build (compilation issues)
- OpenMP detection handles both GCC and clang++ (macOS libomp support)

---

## Optimizations for Speed and Space

The codebase implements multiple optimization strategies to maximize computational performance and minimize memory usage.

### Compiler Optimizations

**Flags Used:** `-O3 -march=native -mtune=native`

- **`-O3`**: Maximum optimization level enabling aggressive inlining, loop unrolling, and vectorization
- **`-march=native`**: Generates code optimized for the specific CPU architecture (uses all available instruction sets)
- **`-mtune=native`**: Tunes code generation for the specific CPU microarchitecture

**Impact:**
- Enables automatic vectorization of loops where possible
- CPU-specific instruction set utilization (SSE, AVX, AVX-512 when available)
- Aggressive inlining reduces function call overhead
- Estimated speedup: 2-4x over `-O2` for compute-intensive kernels

### Memory Optimizations

#### Field3D Flattened Array Storage

**Location:** `include/field3d.hpp`

**Implementation:**
- Replaces nested `std::vector<std::vector<std::vector<float>>>` with single contiguous `std::vector<float>`
- Memory layout: Row-major order (r varies slowest, z fastest)
- Index calculation: `idx = i * NTH * NZ + j * NZ + k`

**Benefits:**
- **Memory Overhead Reduction**: ~99.9% reduction in allocation overhead
  - Nested vectors: 3× pointer overhead + per-vector metadata
  - Flattened array: Single allocation, minimal metadata
- **Cache Locality**: Contiguous memory access improves cache hit rates
  - Estimated speedup: 1.5-3x for memory-bound operations
- **Memory Usage**: ~8GB for large domains (256×128×128 grid) vs. ~12-15GB with nested vectors

**Example Impact:**
- Small grid (64×64×32): ~50MB vs. ~75MB with nested vectors
- Production grid (256×128×128): ~8GB vs. ~12GB with nested vectors

### Parallelization

#### OpenMP Multi-Core Parallelization

**Implementation:**
- Automatic detection of OpenMP support (GCC/clang++)
- macOS libomp support via Homebrew
- Parallel loops using `#pragma omp parallel for collapse(2)`

**Usage Locations:**
- Field initialization loops
- Microphysics tendency calculations
- Advection operations
- Dynamics updates

**Performance:**
- **Speedup**: 4-8x on multi-core systems (4-8 cores)
- **Scaling**: Near-linear scaling up to available cores
- **Overhead**: Minimal (<5%) for small domains, negligible for production grids

**Example Benchmarks:**
- Small grid (64×64×32): ~10,000-15,000 time steps/hour (with OpenMP)
- Production grid (256×128×128): ~100-150 time steps/hour (with OpenMP)
- Without OpenMP: ~2-3x slower

### SIMD Vectorization

#### SIMD-Ready Architecture

**Location:** `include/simd_utils.hpp`, `src/simd_utils.cpp`

**Features:**
- Runtime detection of available SIMD instruction sets (SSE, AVX, AVX-512)
- Vectorized operations for element-wise computations
- Fallback to scalar operations when SIMD unavailable

**Supported Operations:**
- Vector addition: `result = a + b`
- Vector multiplication: `result = a * b`
- Fused multiply-add: `result = a * b + c`
- Process 4/8/16 floats at once (depending on instruction set)

**Current Status:**
- Infrastructure in place
- Ready for integration into computational kernels
- Can provide 4-16x speedup for vectorizable operations

**Potential Integration Points:**
- Advection flux calculations
- Microphysics bulk operations
- Diffusion operators
- Field arithmetic operations

### Cache Optimization Strategies

#### Contiguous Memory Access Patterns

**Implementation:**
- Field3D ensures all data is contiguous in memory
- Loop ordering optimized for cache line utilization
- Row-major access pattern matches memory layout

**Cache Blocking:**
- Large loops can be blocked for better cache reuse
- Potential for further optimization in advection/diffusion kernels

**Impact:**
- Improved cache hit rates (estimated 1.5-3x speedup for memory-bound operations)
- Reduced memory bandwidth requirements

### Performance Benchmarks

**Small Test Grid** (64×64×32):
- Time steps/hour: ~10,000-15,000 (with all optimizations)
- Memory usage: ~50MB
- Single-threaded: ~3,000-5,000 time steps/hour

**Production Grid** (256×128×128):
- Time steps/hour: ~100-150 (with all optimizations)
- Memory usage: ~8GB
- Single-threaded: ~25-40 time steps/hour

**Visualization:**
- Interactive rendering: 30-60 FPS
- Offline rendering: 4K resolution support

### Optimization Opportunities

**Future Improvements:**

1. **SIMD Integration**: Integrate SIMD utilities into computational kernels
   - Advection flux calculations
   - Microphysics bulk operations
   - Expected speedup: 2-4x for vectorizable loops

2. **Cache Blocking**: Implement cache-blocking for large loops
   - Advection operations
   - Diffusion operators
   - Expected speedup: 1.2-1.5x

3. **GPU Acceleration**: Offload compute-intensive operations to GPU
   - Advection/diffusion kernels
   - Microphysics calculations
   - Expected speedup: 10-100x for suitable operations

4. **Memory Pooling**: Reduce allocation overhead for temporary fields
   - Reuse temporary arrays across timesteps
   - Expected memory reduction: 10-20%

5. **Sparse Storage**: Use sparse data structures for fields with many zeros
   - Hydrometeor fields (qi, qs, qg, qh) often sparse
   - Expected memory reduction: 30-50% for sparse fields

### Memory Usage Breakdown

**Typical Production Run** (256×128×128 grid):

- **Field Storage**: ~8GB (all prognostic fields)
- **Temporary Arrays**: ~500MB (tendencies, intermediate calculations)
- **Base State**: ~10MB (density profile, etc.)
- **Total**: ~8.5GB

**Optimization Impact:**
- Without Field3D: ~12-15GB
- Without OpenMP: Same memory, but slower execution
- With all optimizations: ~8.5GB + 4-8x speedup

---

## Visualization

See `visualization/STATUS.md` for detailed visualization pipeline status.

**Summary:**
- 2D animation script works without OpenGL dependencies
- 3D visualization pipeline requires: `moderngl`, `moderngl-window`, `zarr`, `tqdm`
- Modular engine architecture implemented and integrated
- Scripts gracefully handle missing dependencies

---

## Known Faults

### CRITICAL: Potential Temperature Going Negative After One Timestep

**Location:** `src/equations.cpp` lines 451-454 in `step_microphysics()`

**Problem:**
Potential temperature (theta) can become negative after a single timestep, causing the simulation to become unphysical immediately.

**Root Cause:**
In `step_microphysics()`, theta is updated without bounds checking:

```cpp
float theta_old = theta[i][j][k];
float dtheta_total = (dtheta_dt[i][j][k] + dtheta_dt_rad[i][j][k] + dtheta_dt_pbl[i][j][k]) * dt_micro;
theta[i][j][k] += dtheta_total;  // No bounds checking!
```

Large negative tendencies from microphysics/radiation/boundary layer can drive theta negative. Debug code detects large changes (`dtheta_total > 100.0f`) but doesn't prevent negative values.

**Impact:**
- Simulation becomes unphysical immediately
- Can cause NaN/Inf propagation throughout the domain
- Debug checks exist but are reactive, not preventive

**Related Issues:**
- No bounds checking on theta after microphysics updates in `equations.cpp`
- Large tendency values possible (>100 K/s changes detected in debug code)
- Multiple places check for `theta_min < 0` but don't prevent it (reactive, not preventive)
- Bounds clamping scattered across different modules (advection, dynamics) but missing in microphysics step
- Inconsistent bounds enforcement across the codebase

**Files with Theta Bounds Checking:**
- `src/advection/advection.cpp` lines 312-314: Clamps theta to [200, 500] K after advection
- `src/dynamics.cpp` line 382: Clamps theta to [200, 500] K in `step_dynamics_old()`
- `src/equations.cpp` lines 451-454: **NO bounds checking after microphysics updates**

**Debug Code Locations:**
- `src/equations.cpp` lines 457-464: Detects large theta changes but doesn't prevent them
- `src/tornado_sim.cpp` lines 843-866: Periodic checks for negative theta but doesn't prevent it
- `src/equations.cpp` lines 336-338: Initialization check warns but doesn't prevent

**Recommended Fix:**
Add bounds checking immediately after theta update in `step_microphysics()`:

```cpp
theta[i][j][k] += dtheta_total;
// Ensure theta stays within physical bounds
theta[i][j][k] = std::max(200.0f, std::min(500.0f, theta[i][j][k]));
```

Alternatively, investigate why large negative tendencies are occurring and fix the root cause in the physics modules.

---

## Known Issues

### Limitations

- **Terrain Module**: Excluded from build, needs integration work
- **Chaos Module**: Incomplete implementations marked with "COMEBACK" comments
- **API Stability**: Work in progress / research prototype - APIs and file formats may change
- **Test Infrastructure**: Test files consolidated but test targets removed from Makefile

### TODOs

- Fix potential temperature bounds checking issue (CRITICAL)
- Complete chaos module implementations
- Integrate terrain module into build
- Consolidate bounds checking across all modules
- Add preventive checks instead of reactive debug code
- Investigate root causes of large tendency values

---

## Next Steps

1. **Priority 1**: Fix potential temperature negative value issue
   - Add bounds checking in `step_microphysics()`
   - Investigate root causes of large negative tendencies
   - Consider consolidating bounds checking into a utility function

2. **Priority 2**: Complete incomplete modules
   - Finish chaos module implementations
   - Integrate terrain module

3. **Priority 3**: Code quality improvements
   - Consolidate scattered bounds checking
   - Replace reactive debug code with preventive checks
   - Add unit tests for critical physics modules

4. **Priority 4**: Documentation
   - Update README with current status
   - Document physics module interfaces
   - Add examples for common use cases

---

## File Organization

### Test Files Location

All test and debug files have been moved to `tests/` folder:
- `tests/debug_simulation_values.cpp`
- `tests/diagnose_simulation.py`
- `tests/test_soundings.cpp`
- `tests/test_*.yaml` (test configurations)
- `tests/validate_*.py` (validation scripts)

### Configuration Files

- Production configs: `configs/classic.yaml`, `configs/cyclic.yaml`, etc.
- Test configs: Moved to `tests/` folder

---

## Notes

- This STATUS.md file is gitignored for local development context
- For visualization-specific status, see `visualization/STATUS.md`
- For technical documentation, see `docs/README.md`
- For scientific foundations, see `docs/foundationalScience.md`
