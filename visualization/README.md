# SupercellModel 3D Visualization

This folder contains the complete 3D visualization pipeline for supercell thunderstorm simulations.

## Overview

The visualization system provides:
- **Interactive 3D viewer** for real-time exploration of storm structure
- **Offline video renderer** for high-quality animations
- **Data conversion tools** to transform simulation output to visualization format
- **Pipeline runner** for automated workflows

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
   python visualization/render_supercell_3d.py --input data/supercell.zarr --output storm.mp4
   ```

4. **Run complete pipeline:**
   ```bash
   python visualization/run_3d_pipeline.py --all
   ```

## Files

### Core Components
- `supercell_renderer.py` - Interactive 3D OpenGL viewer with volume ray marching
- `render_supercell_3d.py` - Offline video renderer for high-quality animations
- `convert_npy_to_zarr.py` - Data conversion from NPY exports to Zarr format (includes diagnostics)
- `run_3d_pipeline.py` - Automated pipeline runner (includes quick start guide)

### Shaders
- `shaders/volume.vert` - Vertex shader for volume rendering
- `shaders/volume.frag` - Fragment shader with ray marching algorithm

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

### Custom Transfer Functions
Modify `shaders/volume.frag` to customize color mapping and opacity.

### Camera Paths
Extend `supercell_renderer.py` with custom camera trajectories.

### Multi-Field Visualization
Add new fields to the transfer function in the fragment shader.

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
├── qr/        # float32[time, z, y, x] - rain mixing ratio
├── u/, v/, w/ # float32[time, z, y, x] - wind components
├── x, y, z    # float32[x/y/z] - coordinate arrays
└── attrs      # metadata (dx, dy, dz, etc.)
```

## Available Fields

The visualization system can render the following atmospheric fields (when available in your simulation data):

### Primary Fields
- **theta** - Potential temperature (buoyancy/warm updrafts)
- **qv** - Water vapor mixing ratio
- **qc** - Cloud water mixing ratio
- **qr** - Rain water mixing ratio

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
Modify `shaders/volume.frag` to customize how different fields are mapped to colors and opacity.

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

## See Also

- `../docs/foundationalScience.md` - Scientific foundation and literature references
- `../docs/3d_visualization_guide.md` - Complete technical documentation
- `../README.md` - Project overview and setup instructions
- `../configs/` - Example simulation configurations
