"""
Contour rendering pass for 2D/3D contour visualization.
"""

import moderngl as mgl
import numpy as np
from typing import Optional

from ..core.render_engine import RenderPass
from ..core.camera import Camera
from ..core.shader_manager import ShaderManager


class ContourRenderPass(RenderPass):
    """
    Render pass for contour line visualization.
    
    Creates dense contour lines similar to scientific visualization styles.
    """
    
    def __init__(self, shader_manager: ShaderManager,
                 volume_texture: mgl.Texture3D,
                 contour_spacing: float = 0.1,
                 contour_width: float = 0.02,
                 contour_color: tuple = (0.0, 0.0, 0.0),
                 opacity: float = 1.0):
        """
        Initialize contour render pass.
        
        Args:
            shader_manager: Shader manager instance
            volume_texture: 3D texture containing field data
            contour_spacing: Spacing between contour lines
            contour_width: Width of contour lines
            contour_color: RGB color for contours (0-1 range)
            opacity: Opacity of contours
        """
        self.shader_manager = shader_manager
        self.volume_texture = volume_texture
        self.contour_spacing = contour_spacing
        self.contour_width = contour_width
        self.contour_color = np.array(contour_color, dtype=np.float32)
        self.opacity = opacity
        
        self.program: Optional[mgl.Program] = None
        self.cube_vao: Optional[mgl.VertexArray] = None
        self._initialized = False
    
    def _initialize(self, ctx: mgl.Context):
        """Initialize shader program and geometry"""
        if self._initialized:
            return
        
        # Load shader program
        self.program = self.shader_manager.load_shader(
            "contour",
            "contour"
        )
        
        # Create cube VAO (same geometry as volume rendering)
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
    
    def render(self, ctx: mgl.Context, camera: Camera, width: int, height: int, **kwargs):
        """
        Execute contour rendering pass.
        
        Args:
            ctx: ModernGL context
            camera: Camera for view/projection matrices
            width: Viewport width
            height: Viewport height
            **kwargs: Additional arguments
        """
        if not self._initialized:
            self._initialize(ctx)
        
        # Get matrices
        view_matrix = camera.get_view_matrix()
        proj_matrix = camera.get_projection_matrix(width, height)
        mvp_matrix = proj_matrix @ view_matrix
        
        # Set uniforms
        self.program['mvpMatrix'].write(mvp_matrix.tobytes())
        self.program['contourSpacing'].value = self.contour_spacing
        self.program['contourWidth'].value = self.contour_width
        self.program['contourColor'].write(self.contour_color.tobytes())
        self.program['opacity'].value = self.opacity
        
        # Bind volume texture
        self.volume_texture.use(0)
        
        # Render
        self.cube_vao.render(mgl.TRIANGLES)
    
    def set_contour_spacing(self, spacing: float):
        """Set spacing between contour lines"""
        self.contour_spacing = spacing
    
    def set_contour_width(self, width: float):
        """Set width of contour lines"""
        self.contour_width = width
    
    def set_contour_color(self, color: tuple):
        """Set contour color (RGB, 0-1 range)"""
        self.contour_color = np.array(color, dtype=np.float32)
