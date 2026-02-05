# Quick Start Guide - Which Script to Use?

## Simple Decision Tree

```
Do you want to...
├─ Explore data interactively?
│  └─ Use: python visualization/supercell_renderer.py --input data.zarr
│
├─ Create a video animation?
│  ├─ Standard color visualization?
│  │  └─ Use: python visualization/render_supercell_3d.py --input data.zarr --output video.mp4
│  │
│  ├─ Scientific grayscale style?
│  │  └─ Use: python visualization/run_3d_pipeline.py --render-video --scientific --output scientific.mp4
│  │
│  └─ Custom shader?
│     └─ Use: python -m visualization.renderers.offline_renderer --input data.zarr --shader volume_grayscale --output video.mp4
│
├─ Convert simulation data?
│  └─ Use: python visualization/convert_npy_to_zarr.py --input data/exports --output data.zarr
│
└─ Complete workflow (convert + explore + render)?
   └─ Use: python visualization/run_3d_pipeline.py --all
```

## Common Use Cases

### 1. First Time Setup
```bash
# Guided setup for new users
python visualization/run_3d_pipeline.py --quick-start
```

### 2. Interactive Exploration
```bash
# Explore your simulation in 3D
python visualization/supercell_renderer.py --input data/supercell.zarr --field theta

# Try different fields
python visualization/supercell_renderer.py --input data/supercell.zarr --field qr  # Rain
python visualization/supercell_renderer.py --input data/supercell.zarr --field qv  # Water vapor
```

### 3. Create Video Animations

**Standard color visualization:**
```bash
python visualization/render_supercell_3d.py --input data/supercell.zarr --output storm.mp4 --fps 30
```

**Scientific-style (grayscale, side-view):**
```bash
python visualization/run_3d_pipeline.py --render-video --scientific --output scientific.mp4
```

**Custom shader:**
```bash
python -m visualization.renderers.offline_renderer \
    --input data/supercell.zarr \
    --output custom.mp4 \
    --shader volume_grayscale \
    --field theta
```

### 4. Convert Simulation Data
```bash
# Convert NPY exports to Zarr format
python visualization/convert_npy_to_zarr.py \
    --input data/exports \
    --output data/supercell.zarr \
    --max-timesteps 100
```

## Script Reference

| Script | Purpose | When to Use |
|--------|---------|-------------|
| `supercell_renderer.py` | Interactive 3D viewer | Real-time exploration, testing visualizations |
| `render_supercell_3d.py` | Video rendering | Creating standard color animations |
| `run_3d_pipeline.py` | Complete pipeline | Automated workflows, quick-start, scientific style |
| `convert_npy_to_zarr.py` | Data conversion | Converting simulation output to visualization format |
| `create_animation.py` | Simple 2D animation | Quick 2D plots, fallback when OpenGL unavailable |

## Advanced: Using the Modular Engine Directly

For power users who want to customize rendering:

```python
# Use new renderers directly
from visualization.renderers.offline_renderer import OfflineRenderer
from visualization.core.transfer_function import GrayscaleTransferFunction
from visualization.core.camera import ScientificCamera

renderer = OfflineRenderer("data.zarr", shader="volume_grayscale")
renderer.set_transfer_function(GrayscaleTransferFunction("theta"))
renderer.set_camera(ScientificCamera())
renderer.render_animation("output.mp4")
```

See [ENGINE_GUIDE.md](ENGINE_GUIDE.md) for complete documentation.

## Tips

1. **Start with quick-start**: `python visualization/run_3d_pipeline.py --quick-start`
2. **Use pipeline script**: It handles data conversion, viewer, and rendering in one command
3. **Check available fields**: `python visualization/supercell_renderer.py --input data.zarr --list-fields`
4. **For publications**: Use `--scientific` flag for grayscale scientific-style visualizations
5. **Custom shaders**: Place `.frag` files in `shaders/volume/` and reference by name

## Need Help?

- **Quick questions**: See this guide
- **Detailed usage**: See [README.md](README.md)
- **Engine customization**: See [ENGINE_GUIDE.md](ENGINE_GUIDE.md)
- **Scientific background**: See `../docs/foundationalScience.md`
