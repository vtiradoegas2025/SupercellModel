"""
Interactive renderer wrapper for window-based visualization.

This is a simplified version that uses the new engine components.
For full interactive features, the original supercell_renderer.py can be used.
"""

import moderngl_window as mglw
import numpy as np
from pathlib import Path
from typing import Optional

from .base_renderer import BaseRenderer
from ..core.transfer_function import ColorTransferFunction, GrayscaleTransferFunction
from ..core.camera import OrbitCamera
from ..passes.volume_pass import VolumeRenderPass


class InteractiveRenderer(BaseRenderer, mglw.WindowConfig):
    """
    Interactive renderer using moderngl-window for real-time visualization.
    
    Uses the new modular engine components while maintaining window-based interface.
    """
    
    title = "TornadoModel 3D Viewer"
    gl_version = (3, 3)
    aspect_ratio = None
    vsync = True
    
    @classmethod
    def add_arguments(cls, parser):
        parser.add_argument('--input', required=True, help='Data path (Zarr or NPY)')
        parser.add_argument('--field', default='theta', help='Primary field to visualize')
        parser.add_argument('--shader', default='volume_color', 
                          choices=['volume_color', 'volume_grayscale'],
                          help='Shader to use')
    
    def __init__(self, **kwargs):
        # Get arguments from argv if available
        if hasattr(kwargs, 'argv'):
            argv = kwargs['argv']
            data_path = argv.input
            primary_field = getattr(argv, 'field', 'theta')
            shader = getattr(argv, 'shader', 'volume_color')
        else:
            data_path = kwargs.get('input', 'data')
            primary_field = kwargs.get('field', 'theta')
            shader = kwargs.get('shader', 'volume_color')
        
        # Store for later initialization
        self._init_args = {
            'data_path': data_path,
            'primary_field': primary_field,
            'shader': shader
        }
        
        # Initialize window config first
        mglw.WindowConfig.__init__(self, **kwargs)
        
        # Initialize base renderer after window is created
        BaseRenderer.__init__(
            self,
            self._init_args['data_path'],
            self._init_args['primary_field'],
            width=self.window_size[0],
            height=self.window_size[1]
        )
        
        # Setup volume render pass
        self.volume_pass = VolumeRenderPass(
            self.engine.shader_manager,
            self.transfer_function,
            self.volume_texture,
            vertex_shader="volume",
            fragment_shader=self._init_args['shader']
        )
        self.engine.add_render_pass(self.volume_pass)
        
        # Animation state
        self.current_timestep = 0
        self.is_playing = True
        self.last_time = 0.0
        
        # Load initial timestep
        self.load_timestep(0)
    
    def render(self, time: float, frame_time: float):
        """Render frame (called by moderngl-window)"""
        # Update timestep if playing
        if self.is_playing:
            # Simple time-based animation
            timestep = int(time * 0.1) % self.num_timesteps
            if timestep != self.current_timestep:
                self.current_timestep = timestep
                self.load_timestep(timestep)
        
        # Render using engine
        pixels = self.engine.render_frame(timestep=self.current_timestep)
        
        # Display (moderngl-window handles this automatically via FBO)
    
    def key_event(self, key, action, modifiers):
        """Handle keyboard input"""
        if action == self.keys.ACTION_PRESS:
            if key == self.keys.SPACE:
                self.is_playing = not self.is_playing
            elif key == self.keys.R:
                self.current_timestep = 0
                self.load_timestep(0)
            elif key == self.keys.UP:
                self.volume_pass.set_opacity_scale(
                    min(1.0, self.volume_pass.opacity_scale + 0.1)
                )
            elif key == self.keys.DOWN:
                self.volume_pass.set_opacity_scale(
                    max(0.0, self.volume_pass.opacity_scale - 0.1)
                )
    
    def mouse_drag_event(self, x, y, dx, dy):
        """Handle mouse drag for camera rotation"""
        if isinstance(self.camera, OrbitCamera):
            angle_delta = dx * 0.01
            new_angle = self.camera.angle + angle_delta
            self.camera.set_angle(new_angle)
    
    @property
    def num_timesteps(self) -> int:
        """Get number of timesteps"""
        if self.data_manager.format == 'zarr':
            sample_field = self.data_manager.get_all_timesteps(self.primary_field)
            return sample_field.shape[0] if sample_field is not None else 1
        else:
            import glob
            step_dirs = glob.glob(str(self.data_path / "step_*"))
            return len(step_dirs)


def main():
    """Run interactive renderer"""
    mglw.run_window_config(InteractiveRenderer)


if __name__ == "__main__":
    main()
