# TornadoModel 3D Visualization

This folder contains the complete 3D visualization pipeline for supercell thunderstorm simulations.

## Overview

The visualization system provides:
- **Interactive 3D viewer** for real-time exploration of storm structure
- **Offline video renderer** for high-quality animations
- **Data conversion tools** to transform simulation output to visualization format
- **Pipeline runner** for automated workflows
- **Modular rendering engine** with shader support (game-engine-like architecture)

## Architecture

The visualization system is built on a modular, game-engine-like rendering framework:

- **Core Engine** (`core/`): Data management, shader system, rendering engine, cameras
- **Render Passes** (`passes/`): Volume rendering, contour rendering, isosurfaces
- **Renderers** (`renderers/`): High-level renderer implementations
- **Shaders** (`shaders/`): Organized shader library with support for custom shaders

See [ENGINE_GUIDE.md](ENGINE_GUIDE.md) for detailed documentation on the modular system.

## Quick Start

### Fastest Way to Get Started
For new users, run the interactive quick start guide:
```bash
python visualization/run_3d_pipeline.py --quick-start
```
This guided script will check your setup, convert data, and walk you through visualization options.

### Prerequisites
```bash
pip install zarr numpy moderngl moderngl-window pillow tqdm
```

### Manual Setup (Alternative)

1. **Convert simulation data:**
   ```bash
   python visualization/convert_npy_to_zarr.py --input data/exports --output data/supercell.zarr
   ```

2. **Run interactive viewer:**
   ```bash
   python visualization/supercell_renderer.py --input data/supercell.zarr
   ```

3. **Render video:**
   ```bash
   # Standard color visualization
   python visualization/render_supercell_3d.py --input data/supercell.zarr --output storm.mp4
   
   # Scientific-style (grayscale)
   python visualization/run_3d_pipeline.py --render-video --scientific --output scientific.mp4
   
   # Custom shader
   python -m visualization.renderers.offline_renderer --input data/supercell.zarr --shader volume_grayscale --output gray.mp4
   ```

4. **Run complete pipeline:**
   ```bash
   python visualization/run_3d_pipeline.py --all
   ```

## Files

### Entry Points (User-Facing Scripts)
- `supercell_renderer.py` - Interactive 3D OpenGL viewer (uses new engine internally)
- `render_supercell_3d.py` - Offline video renderer (uses new engine internally)
- `run_3d_pipeline.py` - Automated pipeline runner with quick-start guide
- `convert_npy_to_zarr.py` - Data conversion from NPY exports to Zarr format
- `create_animation.py` - Simple 2D matplotlib animation (fallback, no OpenGL required)

### Modular Engine Components
- `core/` - Core engine (DataManager, ShaderManager, RenderEngine, Camera, TransferFunction)
- `renderers/` - High-level renderer implementations
- `passes/` - Render pass implementations (volume, contour)
- `examples/` - Example visualizations (scientific style)

### Shader Library
- `shaders/volume/` - Volume rendering shaders
  - `volume.vert` - Base vertex shader
  - `volume_color.frag` - Color transfer function
  - `volume_grayscale.frag` - Grayscale (scientific style)
- `shaders/contour/` - Contour rendering shaders
- `shaders/common/` - Shared utilities (lighting, utils)

## Controls

### Interactive Viewer
- **Space**: Play/pause animation
- **R**: Reset to beginning
- **W**: Toggle wireframe debug view
- **↑/↓**: Adjust opacity scale
- **Left/Right arrows**: Adjust brightness
- **Mouse drag**: Orbit camera

## Data Format

The system expects simulation data in Zarr format with structure:
```
supercell.zarr/
├── theta/        # Potential temperature (time, z, y, x)
├── qr/           # Rain water mixing ratio (time, z, y, x)
├── u/, v/, w/    # Wind components (time, z, y, x)
├── x, y, z      # Coordinate arrays
└── attrs         # Metadata (dx, dy, dz, etc.)
```

## Pipeline Workflow

1. **Simulation**: Run C++ model with periodic NPY exports
2. **Conversion**: Transform cylindrical data to Cartesian Zarr
3. **Visualization**: Interactive exploration or automated rendering
4. **Export**: High-quality videos for presentations/publications

## Performance

- **GPU Requirements**: OpenGL 3.3+ compatible GPU
- **Memory**: ~4MB per timestep for 128×128×64 volumes
- **Rendering**: 30-60 FPS interactive, 4K offline rendering

## Troubleshooting

### Missing Dependencies
```bash
pip install --upgrade zarr numpy moderngl moderngl-window pillow tqdm
```

### OpenGL Issues
- Ensure graphics drivers are up to date
- For headless rendering: `export DISPLAY=:0` or use EGL backend

### Memory Issues
- Reduce volume resolution in data conversion
- Use downsampling: `--downsample 2 2 2`
- Limit cached timesteps

## Advanced Usage

### Using the Modular Engine

The visualization system now uses a modular engine internally. Advanced users can:

```python
from visualization.renderers.offline_renderer import OfflineRenderer

# Create renderer with custom shader
renderer = OfflineRenderer("data/supercell.zarr", shader="volume_grayscale")
renderer.render_animation("output.mp4")
```

### Custom Shaders

Add custom shaders by placing `.frag` files in `shaders/volume/`:

```bash
# Your custom shader
visualization/shaders/volume/my_custom.frag

# Use it
python -m visualization.renderers.offline_renderer --input data.zarr --shader my_custom
```

Shaders support `#include` directives for shared code:
```glsl
#include "common/utils.glsl"
#include "common/lighting.glsl"
```

### Custom Transfer Functions

Create custom transfer functions:
```python
from visualization.core.transfer_function import CustomTransferFunction

def my_tf(primary_field, secondary_fields, timestep):
    # Your custom logic
    return rgba_data

tf = CustomTransferFunction(my_tf)
renderer.set_transfer_function(tf)
```

### Scientific Style Visualization

For Lewellen et al. (2008) style visualizations:
```bash
python -m visualization.examples.scientific_style --input data.zarr --output scientific.mp4
```

### Camera Paths
Use the Camera system for custom camera trajectories:
```python
from visualization.core.camera import OrbitCamera, ScientificCamera

camera = ScientificCamera()
camera.set_side_view()  # Side view for downdraft visualization
renderer.set_camera(camera)
```

## Troubleshooting

### Zarr Structure Issues

If you get an error like `AttributeError: 'Group' object has no attribute 'shape'`, your zarr file has an incompatible structure. This can happen with older versions of the converter.

**Diagnose the issue:**
```bash
python visualization/convert_npy_to_zarr.py --diagnose --output data/supercell.zarr
```

### Performance Tips

- For large datasets, the renderer loads data on-demand
- Use `--max_timesteps` in the converter to limit data size
- Close the renderer window to free GPU memory

## File Formats

### Input: NPY Exports
The converter expects simulation output in directories like:
```
data/exports/step_000000/
├── th0_theta.npy, th1_theta.npy, ... (128 theta slices)
├── th0_qr.npy, th1_qr.npy, ...      (precipitation)
└── ...
```

### Output: Zarr Arrays
The converter produces Zarr files with structure:
```
storm.zarr/
├── theta/     # float32[time, z, y, x] - potential temperature
├── qv/        # float32[time, z, y, x] - water vapor mixing ratio
├── qc/        # float32[time, z, y, x] - cloud water mixing ratio
├── qr/        # float32[time, z, y, x] - rain mixing ratio
├── qh/        # float32[time, z, y, x] - hail mixing ratio
├── qg/        # float32[time, z, y, x] - graupel mixing ratio
├── u/, v/, w/ # float32[time, z, y, x] - wind components
├── rho/       # float32[time, z, y, x] - density
├── p/         # float32[time, z, y, x] - pressure
├── reflectivity_dbz/ # float32[time, z, y, x] - radar reflectivity (note: exported as "radar")
├── tracer/    # float32[time, z, y, x] - passive tracer
├── x, y, z    # float32[x/y/z] - coordinate arrays
└── attrs      # metadata (dx, dy, dz, etc.)
```

**Note on Field Name Mapping**: The simulation exports radar reflectivity as `radar`, but the converter automatically maps it to `reflectivity_dbz` in the Zarr file for consistency with visualization tools.

## Available Fields

The visualization system can render the following atmospheric fields (when available in your simulation data):

### Primary Fields
- **theta** - Potential temperature (buoyancy/warm updrafts)
- **qv** - Water vapor mixing ratio
- **qc** - Cloud water mixing ratio
- **qr** - Rain water mixing ratio
- **qh** - Hail mixing ratio
- **qg** - Graupel mixing ratio
- **tracer** - Passive tracer field

### Wind Fields
- **u, v, w** - Wind components (radial, azimuthal, vertical)

### Diagnostic Fields (when computed)
- **reflectivity_dbz** - Radar reflectivity
- **buoyancy** - Buoyancy field
- **helicity** - Storm-relative helicity

## Advanced Usage

### Multi-Field Visualization
```bash
# Visualize different primary fields
python supercell_renderer.py --input data/supercell.zarr --field qr
python supercell_renderer.py --input data/supercell.zarr --field qv

# Adjust rendering parameters
python supercell_renderer.py --input data/supercell.zarr --opacity-scale 0.8 --brightness 1.5
```

### Batch Processing
```bash
# Render multiple configurations
python run_3d_pipeline.py --input data/supercell.zarr --output videos/ --formats mp4 webm
```

### Custom Transfer Functions
Modify `shaders/volume/volume_color.frag` or create new shaders in `shaders/volume/` to customize how different fields are mapped to colors and opacity.

See [ENGINE_GUIDE.md](ENGINE_GUIDE.md) for detailed shader customization guide.

## Field Mapping in Shaders

The RGBA channels in the volume texture are mapped as follows:
- **Red (R)**: Primary field (temperature/buoyancy) - warm updrafts = red, cold downdrafts = blue
- **Green (G)**: Precipitation (qr) - white/grey clouds and rain
- **Blue (B)**: Wind speed - cyan for high winds
- **Alpha (A)**: Opacity - based on condensate + minimum visibility

## Troubleshooting

### Missing Fields
If a field you expect isn't available:
1. Check what fields were exported from your simulation
2. Ensure the physics scheme computes that field
3. Use `diagnose_zarr.py` to inspect your data structure

### Performance Issues
- Use `--downsample 2 2 2` in conversion for faster rendering
- Limit `--max_timesteps` to reduce memory usage
- Use lower resolution for interactive exploration

## Field Discovery

The visualization system automatically discovers available fields from simulation output (source-of-truth alignment). No hardcoded field lists - the system adapts to whatever fields your simulation exports.

```python
from visualization.core.data_manager import DataManager

data_mgr = DataManager("data/supercell.zarr")
available_fields = data_mgr.list_available_fields()
print(f"Available: {available_fields}")
```

## See Also

- [ENGINE_GUIDE.md](ENGINE_GUIDE.md) - Complete guide to the modular rendering engine
- [README_ENGINE.md](README_ENGINE.md) - Implementation summary
- `../docs/foundationalScience.md` - Scientific foundation and literature references
- `../docs/3d_visualization_guide.md` - Complete technical documentation
- `../README.md` - Project overview and setup instructions
- `../configs/` - Example simulation configurations
