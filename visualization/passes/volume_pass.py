"""
Volume rendering pass for 3D volume visualization.
"""

import moderngl as mgl
import numpy as np
from typing import Optional, Dict

from ..core.render_engine import RenderPass
from ..core.camera import Camera
from ..core.transfer_function import TransferFunction
from ..core.shader_manager import ShaderManager


class VolumeRenderPass(RenderPass):
    """
    Render pass for volume rendering using ray marching.
    """
    
    def __init__(self, shader_manager: ShaderManager,
                 transfer_function: TransferFunction,
                 volume_texture: mgl.Texture3D,
                 vertex_shader: str = "volume",
                 fragment_shader: str = "volume_color",
                 opacity_scale: float = 0.1,
                 brightness: float = 1.0):
        """
        Initialize volume render pass.
        
        Args:
            shader_manager: Shader manager instance
            transfer_function: Transfer function for field-to-color mapping
            volume_texture: 3D texture containing volume data
            vertex_shader: Name of vertex shader
            fragment_shader: Name of fragment shader
            opacity_scale: Opacity scaling factor
            brightness: Brightness multiplier
        """
        self.shader_manager = shader_manager
        self.transfer_function = transfer_function
        self.volume_texture = volume_texture
        self.vertex_shader = vertex_shader
        self.fragment_shader = fragment_shader
        self.opacity_scale = opacity_scale
        self.brightness = brightness
        
        self.program: Optional[mgl.Program] = None
        self.cube_vao: Optional[mgl.VertexArray] = None
        self._initialized = False
    
    def _initialize(self, ctx: mgl.Context):
        """Initialize shader program and geometry"""
        if self._initialized:
            return
        
        # Load shader program
        self.program = self.shader_manager.load_shader(
            self.vertex_shader,
            self.fragment_shader
        )
        
        # Create cube VAO for volume rendering
        vertices = np.array([
            # Front face
            -1, -1,  1,   1, -1,  1,   1,  1,  1,
            -1, -1,  1,   1,  1,  1,  -1,  1,  1,
            # Right face
             1, -1,  1,   1, -1, -1,   1,  1, -1,
             1, -1,  1,   1,  1, -1,   1,  1,  1,
            # Back face
             1, -1, -1,  -1, -1, -1,  -1,  1, -1,
             1, -1, -1,  -1,  1, -1,   1,  1, -1,
            # Left face
            -1, -1, -1,  -1, -1,  1,  -1,  1,  1,
            -1, -1, -1,  -1,  1,  1,  -1,  1, -1,
            # Top face
            -1,  1,  1,   1,  1,  1,   1,  1, -1,
            -1,  1,  1,   1,  1, -1,  -1,  1, -1,
            # Bottom face
            -1, -1, -1,   1, -1, -1,   1, -1,  1,
            -1, -1, -1,   1, -1,  1,  -1, -1,  1,
        ], dtype=np.float32)
        
        vbo = ctx.buffer(vertices.tobytes())
        self.cube_vao = ctx.vertex_array(
            self.program,
            [(vbo, '3f', 'in_position')]
        )
        
        self._initialized = True
    
    def render(self, ctx: mgl.Context, camera: Camera, width: int, height: int, 
               timestep: int = 0, **kwargs):
        """
        Execute volume rendering pass.
        
        Args:
            ctx: ModernGL context
            camera: Camera for view/projection matrices
            width: Viewport width
            height: Viewport height
            timestep: Timestep index
            **kwargs: Additional arguments (e.g., volume_data, secondary_fields)
        """
        if not self._initialized:
            self._initialize(ctx)
        
        # Get matrices
        view_matrix = camera.get_view_matrix()
        proj_matrix = camera.get_projection_matrix(width, height)
        mvp_matrix = proj_matrix @ view_matrix
        
        # Get camera position
        camera_pos = camera.get_camera_position()
        
        # Set uniforms
        self.program['mvpMatrix'].write(mvp_matrix.tobytes())
        self.program['cameraPos'].write(camera_pos.astype(np.float32).tobytes())
        self.program['opacityScale'].value = self.opacity_scale
        self.program['brightness'].value = self.brightness
        self.program['time'].value = timestep * 0.1
        
        # Bind volume texture
        self.volume_texture.use(0)
        
        # Render
        self.cube_vao.render(mgl.TRIANGLES)
    
    def update_volume_data(self, rgba_data: np.ndarray):
        """
        Update volume texture with new RGBA data.
        
        Args:
            rgba_data: RGBA volume data (z, y, x, 4) as float32
        """
        self.volume_texture.write(rgba_data.tobytes())
    
    def set_opacity_scale(self, opacity_scale: float):
        """Set opacity scaling factor"""
        self.opacity_scale = opacity_scale
    
    def set_brightness(self, brightness: float):
        """Set brightness multiplier"""
        self.brightness = brightness
