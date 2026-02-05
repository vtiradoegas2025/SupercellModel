"""
Base renderer class with common functionality for interactive and offline renderers.
"""

import numpy as np
import moderngl as mgl
from pathlib import Path
from typing import Optional, Dict, List

from ..core.data_manager import DataManager
from ..core.render_engine import RenderEngine
from ..core.transfer_function import TransferFunction, ColorTransferFunction
from ..core.camera import Camera, OrbitCamera
from ..core.shader_manager import ShaderManager
from ..passes.volume_pass import VolumeRenderPass


class BaseRenderer:
    """
    Base class for renderers providing common functionality.
    """
    
    def __init__(self, data_path: str, 
                 primary_field: str = "theta",
                 shader_dir: Optional[Path] = None,
                 width: int = 1920,
                 height: int = 1080):
        """
        Initialize base renderer.
        
        Args:
            data_path: Path to simulation data (Zarr or NPY directory)
            primary_field: Primary field to visualize
            shader_dir: Directory containing shaders
            width: Viewport width
            height: Viewport height
        """
        self.data_path = Path(data_path)
        self.primary_field = primary_field
        self.width = width
        self.height = height
        
        # Initialize data manager
        self.data_manager = DataManager(data_path)
        
        # Initialize render engine
        if shader_dir is None:
            shader_dir = Path(__file__).parent.parent / "shaders"
        self.engine = RenderEngine(width, height, shader_dir=shader_dir)
        self.engine.shader_manager.set_context(self.engine.ctx)
        
        # Initialize transfer function
        self.transfer_function: TransferFunction = ColorTransferFunction(primary_field)
        
        # Volume texture
        self.volume_texture: Optional[mgl.Texture] = None
        self._setup_volume_texture()
        
        # Camera
        self.camera: Camera = OrbitCamera()
        self.engine.set_camera(self.camera)
        
        # Render pass
        self.volume_pass: Optional[VolumeRenderPass] = None
    
    def _setup_volume_texture(self):
        """Setup 3D volume texture"""
        nx, ny, nz = self.data_manager.get_grid_dims()
        
        self.volume_texture = self.engine.ctx.texture3d(
            (nx, ny, nz), 4, dtype='f4'
        )
        self.volume_texture.filter = (mgl.LINEAR, mgl.LINEAR)
    
    def load_timestep(self, timestep: int = 0):
        """Load data for a specific timestep and update volume texture"""
        # Get primary field
        primary_data = self.data_manager.get_field(self.primary_field, timestep)
        
        # Get secondary fields
        secondary_fields = {}
        available_fields = self.data_manager.list_available_fields()
        for field_name in available_fields:
            if field_name != self.primary_field:
                try:
                    field_data = self.data_manager.get_field(field_name, timestep)
                    if field_data is not None:
                        secondary_fields[field_name] = field_data
                except:
                    pass
        
        # Apply transfer function
        rgba_data = self.transfer_function.apply(
            primary_data,
            secondary_fields,
            timestep
        )
        
        # Update volume texture
        if self.volume_texture is not None:
            self.volume_texture.write(rgba_data.tobytes())
    
    def set_transfer_function(self, transfer_function: TransferFunction):
        """Set transfer function"""
        self.transfer_function = transfer_function
    
    def set_camera(self, camera: Camera):
        """Set camera"""
        self.camera = camera
        self.engine.set_camera(camera)
    
    def get_available_fields(self) -> List[str]:
        """Get list of available fields"""
        return self.data_manager.list_available_fields()
    
    def get_grid_dims(self) -> tuple:
        """Get grid dimensions"""
        return self.data_manager.get_grid_dims()
