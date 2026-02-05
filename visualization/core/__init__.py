"""
Core rendering engine components for TornadoModel visualization.

This module provides the foundational components for a modular, game-engine-like
rendering system with shader support.
"""

from .data_manager import DataManager
from .shader_manager import ShaderManager
from .transfer_function import (
    TransferFunction,
    ColorTransferFunction,
    GrayscaleTransferFunction,
    CustomTransferFunction
)
from .render_engine import RenderEngine
from .camera import Camera, OrbitCamera, ScientificCamera

__all__ = [
    'DataManager',
    'ShaderManager',
    'TransferFunction',
    'ColorTransferFunction',
    'GrayscaleTransferFunction',
    'CustomTransferFunction',
    'RenderEngine',
    'Camera',
    'OrbitCamera',
    'ScientificCamera',
]
