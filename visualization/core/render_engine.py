"""
Core rendering engine for OpenGL-based volume rendering.

Manages OpenGL context, framebuffers, and rendering pipeline.
"""

import moderngl as mgl
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Any, Union
from abc import ABC, abstractmethod

from .camera import Camera
from .shader_manager import ShaderManager


class RenderPass(ABC):
    """Abstract base class for render passes"""
    
    @abstractmethod
    def render(self, ctx: mgl.Context, camera: Camera, width: int, height: int, **kwargs):
        """Execute render pass"""
        pass


class RenderEngine:
    """
    Core rendering engine managing OpenGL context and rendering pipeline.
    """
    
    def __init__(self, width: int = 1920, height: int = 1080, 
                 ctx: Optional[mgl.Context] = None,
                 shader_dir: Optional[Path] = None):
        """
        Initialize render engine.
        
        Args:
            width: Viewport width
            height: Viewport height
            ctx: ModernGL context (if None, creates standalone context)
            shader_dir: Directory containing shaders
        """
        self.width = width
        self.height = height
        
        # Create or use provided context
        if ctx is None:
            try:
                # Try standard standalone context
                self.ctx = mgl.create_standalone_context()
            except Exception as e1:
                try:
                    # Try with explicit version
                    self.ctx = mgl.create_standalone_context(require=330)
                except Exception as e2:
                    # macOS limitation: requires display for OpenGL context
                    # The interactive viewer works because it creates a visible window
                    raise RuntimeError(
                        f"macOS cannot create headless OpenGL context. "
                        f"Use the interactive viewer: "
                        f"python visualization/supercell_renderer.py --input <data>"
                    )
        else:
            self.ctx = ctx
        
        # Enable OpenGL features
        self.ctx.enable(mgl.DEPTH_TEST)
        self.ctx.enable(mgl.BLEND)
        self.ctx.blend_func = mgl.SRC_ALPHA, mgl.ONE_MINUS_SRC_ALPHA
        
        # Initialize shader manager
        if shader_dir is None:
            shader_dir = Path(__file__).parent.parent / "shaders"
        self.shader_manager = ShaderManager(shader_dir, self.ctx)
        
        # Create framebuffer for offscreen rendering
        self.fbo = self.ctx.framebuffer(
            color_attachments=[self.ctx.texture((width, height), 4)],
            depth_attachment=self.ctx.depth_texture((width, height))
        )
        
        # Render passes
        self.render_passes: List[RenderPass] = []
        
        # Camera
        self.camera: Optional[Camera] = None
    
    def set_camera(self, camera: Camera):
        """Set active camera"""
        self.camera = camera
    
    def add_render_pass(self, render_pass: RenderPass):
        """Add a render pass to the pipeline"""
        self.render_passes.append(render_pass)
    
    def clear_render_passes(self):
        """Clear all render passes"""
        self.render_passes.clear()
    
    def render_frame(self, **kwargs) -> np.ndarray:
        """
        Render a single frame.
        
        Args:
            **kwargs: Additional arguments passed to render passes
            
        Returns:
            Rendered image as numpy array (height, width, 4) uint8
        """
        if self.camera is None:
            raise RuntimeError("Camera not set. Call set_camera() first.")
        
        # Bind framebuffer
        self.fbo.use()
        
        # Clear
        self.ctx.clear(0.1, 0.1, 0.2, 1.0)
        
        # Execute render passes
        for render_pass in self.render_passes:
            render_pass.render(
                self.ctx,
                self.camera,
                self.width,
                self.height,
                **kwargs
            )
        
        # Read back pixels
        pixels = self.fbo.read(components=4, dtype='f4')
        pixels = pixels.reshape((self.height, self.width, 4))
        
        # Convert to 8-bit RGBA
        pixels = np.clip(pixels * 255, 0, 255).astype(np.uint8)
        
        return pixels
    
    def resize(self, width: int, height: int):
        """Resize viewport"""
        self.width = width
        self.height = height
        
        # Recreate framebuffer
        self.fbo.release()
        self.fbo = self.ctx.framebuffer(
            color_attachments=[self.ctx.texture((width, height), 4)],
            depth_attachment=self.ctx.depth_texture((width, height))
        )
    
    def cleanup(self):
        """Clean up resources"""
        if hasattr(self, 'fbo'):
            self.fbo.release()
        self.shader_manager.clear_cache()
