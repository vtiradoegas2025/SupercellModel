"""
Data Manager for loading simulation data.

Source of truth: The simulation code (src/tornado_sim.cpp) exports fields via write_all_fields().
This manager automatically discovers available fields from the actual data files rather than
hardcoding field lists, ensuring visualization stays in sync with simulation.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import glob
import warnings

try:
    import zarr
    ZARR_AVAILABLE = True
except ImportError:
    ZARR_AVAILABLE = False


class DataManager:
    """
    Unified interface for loading simulation data from Zarr or NPY formats.
    
    Automatically discovers available fields from data files, ensuring alignment
    with what the simulation actually exports.
    """
    
    # Fields exported by simulation (from src/tornado_sim.cpp write_all_fields)
    # These are the standard fields, but we auto-detect from actual data
    KNOWN_FIELDS = {
        'u', 'v', 'w',  # Wind components
        'rho', 'p',     # Density, pressure
        'theta',        # Potential temperature
        'qv', 'qc', 'qr', 'qh', 'qg', 'qi', 'qs',  # Moisture fields
        'radar', 'radar_reflectivity',  # Radar (exported as "radar", may be mapped to "radar_reflectivity")
        'tracer',       # Passive tracer
    }
    
    def __init__(self, data_path: Union[str, Path], format: Optional[str] = None):
        """
        Initialize DataManager.
        
        Args:
            data_path: Path to Zarr file or directory containing NPY files
            format: 'zarr' or 'npy'. If None, auto-detect from path
        """
        self.data_path = Path(data_path)
        
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data path does not exist: {self.data_path}")
        
        # Auto-detect format
        if format is None:
            # Zarr stores can be directories or files
            # First check if it's a Zarr store (by trying to open it)
            if ZARR_AVAILABLE:
                try:
                    store = zarr.open(str(self.data_path), mode='r')
                    # Check if it has keys (store.keys() is a generator, convert to list)
                    keys = list(store.keys())
                    if len(keys) > 0:
                        format = 'zarr'
                    else:
                        # Empty zarr, check for NPY structure
                        if self.data_path.is_dir() and any(self.data_path.glob('step_*')):
                            format = 'npy'
                        else:
                            format = 'zarr'  # Default to zarr even if empty
                except:
                    # Not a zarr store, check for NPY structure
                    if self.data_path.is_dir() and any(self.data_path.glob('step_*')):
                        format = 'npy'
                    elif self.data_path.is_file() and self.data_path.suffix == '.zarr':
                        format = 'zarr'
                    else:
                        format = 'npy'  # Default fallback
            else:
                # Zarr not available, check for NPY structure
                if self.data_path.is_dir() and any(self.data_path.glob('step_*')):
                    format = 'npy'
                else:
                    format = 'npy'  # Default fallback
        
        self.format = format
        self.store = None
        self.npy_base_dir = None
        
        # Field cache
        self._field_cache: Dict[str, np.ndarray] = {}
        self._discovered_fields: Optional[List[str]] = None
        self._grid_dims: Optional[Tuple[int, int, int]] = None
        
        # Load data source
        if self.format == 'zarr':
            self._load_zarr()
        else:
            self._load_npy()
        
        # Discover available fields
        self.discover_fields()
    
    def _load_zarr(self):
        """Load Zarr store"""
        if not ZARR_AVAILABLE:
            raise RuntimeError("zarr package not available. Install with: pip install zarr")
        
        self.store = zarr.open(str(self.data_path), mode='r')
        
        # Validate structure
        if 'theta' not in self.store and len(self.store.keys()) == 0:
            raise ValueError(f"Zarr file appears empty: {self.data_path}")
    
    def _load_npy(self):
        """Load NPY directory structure"""
        if not self.data_path.is_dir():
            raise ValueError(f"NPY format requires a directory: {self.data_path}")
        
        # Find step directories
        step_dirs = sorted(glob.glob(str(self.data_path / "step_*")))
        if not step_dirs:
            raise ValueError(f"No step_* directories found in {self.data_path}")
        
        self.npy_base_dir = self.data_path
        self._step_dirs = step_dirs
    
    def discover_fields(self) -> List[str]:
        """
        Automatically discover available fields from data files.
        
        This is the source-of-truth discovery - we find what fields actually exist
        in the data, not what we expect.
        
        Returns:
            List of discovered field names
        """
        if self._discovered_fields is not None:
            return self._discovered_fields
        
        discovered = set()
        
        if self.format == 'zarr':
            # Discover from Zarr keys
            for key in self.store.keys():
                # Skip coordinate arrays and metadata
                if key in ['x', 'y', 'z', 'attrs']:
                    continue
                # Check if it's an array
                item = self.store[key]
                if hasattr(item, 'shape') and len(item.shape) == 4:  # (time, z, y, x)
                    discovered.add(key)
        else:
            # Discover from NPY files in first timestep
            if self._step_dirs:
                first_step = Path(self._step_dirs[0])
                # Look for files matching pattern th*_*.npy
                for npy_file in first_step.glob("th*_*.npy"):
                    # Extract field name: th0_theta.npy -> theta
                    parts = npy_file.stem.split("_", 1)
                    if len(parts) == 2:
                        field_name = parts[1]
                        discovered.add(field_name)
        
        self._discovered_fields = sorted(list(discovered))
        return self._discovered_fields
    
    def list_available_fields(self) -> List[str]:
        """Get list of all available fields"""
        return self.discover_fields()
    
    def has_field(self, field_name: str) -> bool:
        """Check if a field is available"""
        return field_name in self.discover_fields()
    
    def get_grid_dims(self) -> Tuple[int, int, int]:
        """
        Get grid dimensions (nx, ny, nz).
        
        Returns:
            Tuple of (nx, ny, nz)
        """
        if self._grid_dims is not None:
            return self._grid_dims
        
        # Try to get from a known field
        if self.has_field('theta'):
            sample = self.get_field('theta', timestep=0)
            if sample is not None:
                # Shape is (z, y, x) for a single timestep
                self._grid_dims = (sample.shape[2], sample.shape[1], sample.shape[0])
                return self._grid_dims
        
        # Fallback: try to infer from coordinates or metadata
        if self.format == 'zarr' and 'x' in self.store:
            nx = len(self.store['x'])
            ny = len(self.store['y']) if 'y' in self.store else nx
            nz = len(self.store['z']) if 'z' in self.store else 128
            self._grid_dims = (nx, ny, nz)
            return self._grid_dims
        
        # Default fallback
        warnings.warn("Could not determine grid dimensions, using defaults")
        self._grid_dims = (128, 128, 128)
        return self._grid_dims
    
    def get_field(self, field_name: str, timestep: int = 0) -> Optional[np.ndarray]:
        """
        Get field data for a specific timestep.
        
        Args:
            field_name: Name of the field
            timestep: Timestep index
            
        Returns:
            3D array (z, y, x) or None if field not found
        """
        # Check cache first
        cache_key = f"{field_name}_{timestep}"
        if cache_key in self._field_cache:
            return self._field_cache[cache_key]
        
        # Validate field exists
        if not self.has_field(field_name):
            available = self.list_available_fields()
            raise ValueError(
                f"Field '{field_name}' not found. Available fields: {available}"
            )
        
        data = None
        
        if self.format == 'zarr':
            # Load from Zarr
            field_data = self.store[field_name]
            if hasattr(field_data, 'shape'):
                # Direct array access
                if timestep < field_data.shape[0]:
                    data = field_data[timestep]  # (z, y, x)
                else:
                    raise IndexError(f"Timestep {timestep} out of range (max: {field_data.shape[0]-1})")
        else:
            # Load from NPY files
            data = self._load_npy_field(field_name, timestep)
        
        # Cache result
        if data is not None:
            self._field_cache[cache_key] = data
        
        return data
    
    def _load_npy_field(self, field_name: str, timestep: int) -> np.ndarray:
        """Load field from NPY files for a specific timestep"""
        if timestep >= len(self._step_dirs):
            raise IndexError(f"Timestep {timestep} out of range (max: {len(self._step_dirs)-1})")
        
        step_dir = Path(self._step_dirs[timestep])
        
        # Find all theta slices for this field
        slices = []
        theta_idx = 0
        
        while True:
            npy_file = step_dir / f"th{theta_idx}_{field_name}.npy"
            if not npy_file.exists():
                break
            
            slice_data = np.load(npy_file)
            slices.append(slice_data)
            theta_idx += 1
        
        if not slices:
            return None
        
        # Stack slices: each slice is (NZ, NR), stacking gives (NZ, NTH, NR)
        # Then transpose to (NR, NTH, NZ) and convert to Cartesian (NX, NY, NZ)
        stacked = np.stack(slices, axis=1)  # (NZ, NTH, NR)
        data_3d = np.transpose(stacked, (2, 1, 0))  # (NR, NTH, NZ)
        
        # For now, return cylindrical data - conversion to Cartesian happens in renderer
        # TODO: Add cylindrical-to-Cartesian conversion option
        return data_3d
    
    def get_all_timesteps(self, field_name: str) -> Optional[np.ndarray]:
        """
        Get all timesteps for a field.
        
        Args:
            field_name: Name of the field
            
        Returns:
            4D array (time, z, y, x) or None if field not found
        """
        if not self.has_field(field_name):
            return None
        
        if self.format == 'zarr':
            field_data = self.store[field_name]
            if hasattr(field_data, 'shape'):
                return field_data[:]  # (time, z, y, x)
        else:
            # Load all timesteps from NPY
            all_data = []
            for timestep in range(len(self._step_dirs)):
                data = self.get_field(field_name, timestep)
                if data is not None:
                    all_data.append(data)
            
            if all_data:
                return np.stack(all_data, axis=0)  # (time, z, y, x)
        
        return None
    
    def get_coordinates(self) -> Dict[str, np.ndarray]:
        """
        Get coordinate arrays.
        
        Returns:
            Dictionary with 'x', 'y', 'z' arrays
        """
        coords = {}
        
        if self.format == 'zarr':
            for axis in ['x', 'y', 'z']:
                if axis in self.store:
                    coords[axis] = np.array(self.store[axis])
        else:
            # For NPY, we need to infer from grid dimensions
            nx, ny, nz = self.get_grid_dims()
            # Default: assume 1000m resolution, 100m vertical
            coords['x'] = np.linspace(-nx*500, nx*500, nx)
            coords['y'] = np.linspace(-ny*500, ny*500, ny)
            coords['z'] = np.linspace(0, nz*100, nz)
        
        return coords
    
    def clear_cache(self):
        """Clear field cache"""
        self._field_cache.clear()
    
    def get_metadata(self) -> Dict:
        """Get metadata about the dataset"""
        metadata = {
            'format': self.format,
            'path': str(self.data_path),
            'fields': self.list_available_fields(),
            'grid_dims': self.get_grid_dims(),
        }
        
        if self.format == 'zarr' and 'attrs' in self.store:
            metadata['attrs'] = dict(self.store['attrs'])
        
        return metadata
