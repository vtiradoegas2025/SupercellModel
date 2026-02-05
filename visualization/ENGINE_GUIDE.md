# Rendering Engine User Guide

## Overview

The TornadoModel visualization system is built on a modular, game-engine-like rendering framework that supports:

- **Simple shader swapping**: Drop shader files into the `shaders/` directory
- **Advanced shader composition**: Combine multiple shader effects using `#include` directives
- **Multiple rendering modes**: Volume, contour, isosurface rendering
- **Extensible architecture**: Easy to add custom components

## Architecture

### Core Components

1. **DataManager**: Loads simulation data (Zarr/NPY) and auto-discovers available fields
2. **ShaderManager**: Loads, compiles, and composes GLSL shaders
3. **TransferFunction**: Converts simulation fields to RGBA volume textures
4. **RenderEngine**: Manages OpenGL context and rendering pipeline
5. **Camera**: Provides different viewing angles (orbit, scientific, etc.)
6. **RenderPasses**: Modular rendering techniques (volume, contour, etc.)

## Quick Start

### Basic Usage

```python
from visualization.renderers.offline_renderer import OfflineRenderer

# Create renderer
renderer = OfflineRenderer(
    data_path="data/simulation.zarr",
    primary_field="theta",
    shader="volume_color"
)

# Render animation
renderer.render_animation("output.mp4", fps=30)
```

### Command Line

```bash
# Render with color shader
python -m visualization.renderers.offline_renderer \
    --input data/simulation.zarr \
    --output animation.mp4 \
    --shader volume_color

# Render with grayscale (scientific style)
python -m visualization.renderers.offline_renderer \
    --input data/simulation.zarr \
    --output scientific.mp4 \
    --shader volume_grayscale \
    --camera scientific
```

## Adding Custom Shaders

### Simple Shader Swapping

1. Create a new fragment shader in `visualization/shaders/volume/`:

```glsl
// my_custom_shader.frag
#version 330

in vec3 worldPos;
in vec3 viewDir;
out vec4 fragColor;

uniform sampler3D volumeTexture;
uniform vec3 cameraPos;
uniform float opacityScale;
uniform float brightness;

vec4 transferFunction(vec4 sample) {
    // Your custom transfer function
    float intensity = sample.r;
    vec3 color = vec3(intensity, 0.0, 1.0 - intensity);  // Purple gradient
    return vec4(color, intensity * opacityScale);
}

void main() {
    // Standard ray marching code (copy from volume_color.frag)
    // ...
}
```

2. Use it:

```python
renderer = OfflineRenderer(data_path, shader="my_custom_shader")
```

### Using Shader Includes

Share common code between shaders:

```glsl
// my_shader.frag
#include "common/utils.glsl"
#include "common/lighting.glsl"

vec4 transferFunction(vec4 sample) {
    float intensity = rgbToGrayscale(sample.rgb);
    // Use shared utilities
    return vec4(intensity, intensity, intensity, sample.a);
}
```

## Transfer Functions

### Built-in Transfer Functions

- **ColorTransferFunction**: Multi-field color visualization
- **GrayscaleTransferFunction**: Scientific-style grayscale

### Custom Transfer Function

```python
from visualization.core.transfer_function import TransferFunction, CustomTransferFunction

def my_transfer_func(primary_field, secondary_fields, timestep):
    # Your custom logic
    intensity = np.clip(primary_field / 300.0, 0.0, 1.0)
    rgba = np.stack([
        intensity,
        intensity * 0.5,
        1.0 - intensity,
        intensity * 0.8
    ], axis=-1)
    return rgba.astype(np.float32)

tf = CustomTransferFunction(my_transfer_func)
renderer.set_transfer_function(tf)
```

## Camera System

### Orbit Camera

```python
from visualization.core.camera import OrbitCamera

camera = OrbitCamera(
    distance=50.0,
    height=20.0,
    angle=0.0  # radians
)
renderer.set_camera(camera)
```

### Scientific Camera

```python
from visualization.core.camera import ScientificCamera

camera = ScientificCamera()
camera.set_side_view()  # Side view for downdraft visualization
renderer.set_camera(camera)
```

## Render Passes

### Volume Rendering

```python
from visualization.passes.volume_pass import VolumeRenderPass

volume_pass = VolumeRenderPass(
    shader_manager,
    transfer_function,
    volume_texture,
    fragment_shader="volume_grayscale"
)
engine.add_render_pass(volume_pass)
```

### Contour Rendering

```python
from visualization.passes.contour_pass import ContourRenderPass

contour_pass = ContourRenderPass(
    shader_manager,
    volume_texture,
    contour_spacing=0.05,
    contour_width=0.01,
    contour_color=(0.0, 0.0, 0.0)  # Black
)
engine.add_render_pass(contour_pass)
```

## Field Discovery

The DataManager automatically discovers available fields from simulation output:

```python
from visualization.core.data_manager import DataManager

data_mgr = DataManager("data/simulation.zarr")
available_fields = data_mgr.list_available_fields()
print(f"Available fields: {available_fields}")

# Fields are discovered from actual simulation exports
# No hardcoded field lists - stays in sync with simulation
```

## Examples

### Scientific Style Visualization

```python
from visualization.examples.scientific_style import ScientificStyleRenderer

renderer = ScientificStyleRenderer("data/simulation.zarr")
renderer.render_animation("scientific.mp4", fps=10)
```

## Advanced: Shader Composition

For power users, compose multiple shader effects:

```python
from visualization.core.shader_manager import ShaderManager

shader_mgr = ShaderManager("shaders/")

# Compose multiple shader parts
composed = shader_mgr.compose_shaders([
    "volume_base",
    "lighting_effects",
    "post_process"
], shader_type='frag')
```

## Performance Tips

1. **Pre-compute transfer functions**: Cache RGBA volumes for repeated rendering
2. **Use appropriate resolution**: Lower resolution for faster iteration
3. **Limit timesteps**: Use `--max-timesteps` for testing
4. **GPU memory**: Large volumes may require texture size limits

## Troubleshooting

### Shader Not Found

- Check shader file exists in `shaders/volume/` or `shaders/contour/`
- Verify shader name matches filename (without extension)
- Check shader syntax with GLSL validator

### Field Not Found

- Use `data_mgr.list_available_fields()` to see what's actually available
- Check simulation actually exports the field you're requesting
- Field names are case-sensitive

### OpenGL Errors

- Ensure graphics drivers are up to date
- Check OpenGL version (requires 3.3+)
- For headless rendering, may need EGL backend

## API Reference

See inline documentation in:
- `visualization/core/data_manager.py`
- `visualization/core/shader_manager.py`
- `visualization/core/transfer_function.py`
- `visualization/core/render_engine.py`
- `visualization/core/camera.py`
