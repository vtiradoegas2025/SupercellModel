"""
Transfer Function system for converting simulation fields to RGBA volume textures.

Transfer functions map scalar field values to colors and opacity for volume rendering.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple


class TransferFunction(ABC):
    """
    Abstract base class for transfer functions.
    
    Transfer functions convert simulation field data into RGBA volume textures
    suitable for volume rendering.
    """
    
    @abstractmethod
    def apply(self, primary_field: np.ndarray, 
              secondary_fields: Optional[Dict[str, np.ndarray]] = None,
              timestep: int = 0) -> np.ndarray:
        """
        Apply transfer function to field data.
        
        Args:
            primary_field: Primary field data (3D array: z, y, x)
            secondary_fields: Optional dictionary of secondary fields
            timestep: Timestep index (for time-varying transfer functions)
            
        Returns:
            RGBA volume data (z, y, x, 4) as float32
        """
        pass
    
    def normalize_field(self, field: np.ndarray, 
                       vmin: Optional[float] = None,
                       vmax: Optional[float] = None) -> np.ndarray:
        """
        Normalize field to [0, 1] range.
        
        Args:
            field: Input field
            vmin: Minimum value (if None, use field min)
            vmax: Maximum value (if None, use field max)
            
        Returns:
            Normalized field
        """
        if vmin is None:
            vmin = np.nanmin(field)
        if vmax is None:
            vmax = np.nanmax(field)
        
        if vmax > vmin:
            return np.clip((field - vmin) / (vmax - vmin), 0.0, 1.0)
        else:
            return np.zeros_like(field)


class ColorTransferFunction(TransferFunction):
    """
    Color transfer function for multi-field visualization.
    
    Maps:
    - Primary field (theta/temperature) -> Red channel (warm=red, cold=blue)
    - Condensate fields (qr, qc, etc.) -> Green channel (white/grey)
    - Water vapor/wind -> Blue channel (cyan/magenta)
    - Opacity based on condensate and minimum visibility
    """
    
    def __init__(self, primary_field: str = "theta",
                 secondary_fields: Optional[List[str]] = None):
        """
        Initialize color transfer function.
        
        Args:
            primary_field: Name of primary field to visualize
            secondary_fields: List of secondary field names to include
        """
        self.primary_field = primary_field
        self.secondary_fields = secondary_fields or []
    
    def apply(self, primary_field: np.ndarray,
              secondary_fields: Optional[Dict[str, np.ndarray]] = None,
              timestep: int = 0) -> np.ndarray:
        """Apply color transfer function"""
        # Normalize primary field
        primary_min, primary_max = np.nanmin(primary_field), np.nanmax(primary_field)
        if primary_max > primary_min:
            primary_norm = (primary_field - primary_min) / (primary_max - primary_min)
        else:
            primary_norm = np.zeros_like(primary_field)
        
        # Initialize RGBA
        r = primary_norm.copy()
        g = np.zeros_like(primary_norm)
        b = np.zeros_like(primary_norm)
        a = np.zeros_like(primary_norm)
        
        if secondary_fields is None:
            secondary_fields = {}
        
        # Add condensate fields (clouds, rain, ice, etc.)
        condensate_opacity = np.zeros_like(a)
        
        # Rain water mixing ratio
        if 'qr' in secondary_fields:
            qr = secondary_fields['qr']
            qr_norm = np.clip(qr * 1000.0, 0.0, 1.0)
            g += qr_norm * 0.6
            condensate_opacity += qr_norm * 0.5
        
        # Cloud water mixing ratio
        if 'qc' in secondary_fields:
            qc = secondary_fields['qc']
            qc_norm = np.clip(qc * 1000.0, 0.0, 1.0)
            g += qc_norm * 0.4
            condensate_opacity += qc_norm * 0.3
        
        # Cloud ice mixing ratio
        if 'qi' in secondary_fields:
            qi = secondary_fields['qi']
            qi_norm = np.clip(qi * 1000.0, 0.0, 1.0)
            g += qi_norm * 0.3
            condensate_opacity += qi_norm * 0.4
        
        # Snow mixing ratio
        if 'qs' in secondary_fields:
            qs = secondary_fields['qs']
            qs_norm = np.clip(qs * 1000.0, 0.0, 1.0)
            g += qs_norm * 0.2
            condensate_opacity += qs_norm * 0.3
        
        # Hail/graupel mixing ratio
        if 'qh' in secondary_fields:
            qh = secondary_fields['qh']
            qh_norm = np.clip(qh * 1000.0, 0.0, 1.0)
            g += qh_norm * 0.5
            condensate_opacity += qh_norm * 0.6
        
        if 'qg' in secondary_fields:
            qg = secondary_fields['qg']
            qg_norm = np.clip(qg * 1000.0, 0.0, 1.0)
            g += qg_norm * 0.4
            condensate_opacity += qg_norm * 0.5
        
        # Water vapor field
        if 'qv' in secondary_fields:
            qv = secondary_fields['qv']
            qv_norm = np.clip((qv - 0.005) / 0.015, 0.0, 1.0)
            b += qv_norm * 0.5
        
        # Wind field visualization
        if all(k in secondary_fields for k in ['u', 'v', 'w']):
            u = secondary_fields['u']
            v = secondary_fields['v']
            w = secondary_fields['w']
            
            wind_speed = np.sqrt(u**2 + v**2 + w**2)
            wind_norm = np.clip(wind_speed / 50.0, 0.0, 1.0)
            
            b = np.maximum(b, wind_norm)
            condensate_opacity = np.maximum(condensate_opacity, wind_norm * 0.2)
        
        # Add condensate opacity to final alpha
        a = np.maximum(a, condensate_opacity)
        
        # Minimum opacity for visibility
        a = np.maximum(a, 0.05)
        
        return np.stack([r, g, b, a], axis=-1).astype(np.float32)


class GrayscaleTransferFunction(TransferFunction):
    """
    Grayscale transfer function for scientific-style visualization.
    
    Maps scalar field values to grayscale intensity, suitable for
    scientific publications and presentations.
    """
    
    def __init__(self, field_name: str = "theta",
                 invert: bool = False,
                 contrast: float = 1.0):
        """
        Initialize grayscale transfer function.
        
        Args:
            field_name: Name of field to visualize
            invert: If True, invert the grayscale (dark=high, light=low)
            contrast: Contrast multiplier (>1 = higher contrast)
        """
        self.field_name = field_name
        self.invert = invert
        self.contrast = contrast
    
    def apply(self, primary_field: np.ndarray,
              secondary_fields: Optional[Dict[str, np.ndarray]] = None,
              timestep: int = 0) -> np.ndarray:
        """Apply grayscale transfer function"""
        # Normalize field
        field_min, field_max = np.nanmin(primary_field), np.nanmax(primary_field)
        if field_max > field_min:
            intensity = (primary_field - field_min) / (field_max - field_min)
        else:
            intensity = np.zeros_like(primary_field)
        
        # Apply contrast
        intensity = np.power(intensity, 1.0 / self.contrast)
        
        # Invert if requested
        if self.invert:
            intensity = 1.0 - intensity
        
        # Create grayscale RGBA
        r = intensity
        g = intensity
        b = intensity
        
        # Opacity based on intensity (higher values more opaque)
        # For scientific visualization, we want smooth opacity
        a = intensity * 0.8 + 0.1  # Range: 0.1 to 0.9
        
        return np.stack([r, g, b, a], axis=-1).astype(np.float32)


class CustomTransferFunction(TransferFunction):
    """
    User-extensible transfer function.
    
    Users can subclass this or provide a custom function.
    """
    
    def __init__(self, transfer_func):
        """
        Initialize custom transfer function.
        
        Args:
            transfer_func: Callable that takes (primary_field, secondary_fields, timestep)
                          and returns RGBA array
        """
        self.transfer_func = transfer_func
    
    def apply(self, primary_field: np.ndarray,
              secondary_fields: Optional[Dict[str, np.ndarray]] = None,
              timestep: int = 0) -> np.ndarray:
        """Apply custom transfer function"""
        return self.transfer_func(primary_field, secondary_fields or {}, timestep)
