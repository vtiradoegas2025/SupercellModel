"""
Camera system for 3D visualization.

Provides different camera types for various viewing angles and animations.
"""

import numpy as np
import math
from abc import ABC, abstractmethod
from typing import Optional, Tuple, List


class Camera(ABC):
    """Abstract base class for cameras"""
    
    @abstractmethod
    def get_view_matrix(self) -> np.ndarray:
        """Get view matrix (4x4)"""
        pass
    
    @abstractmethod
    def get_projection_matrix(self, width: int, height: int) -> np.ndarray:
        """Get projection matrix (4x4)"""
        pass
    
    @abstractmethod
    def get_camera_position(self) -> np.ndarray:
        """Get camera position in world space"""
        pass
    
    @staticmethod
    def look_at(eye: np.ndarray, center: np.ndarray, up: np.ndarray) -> np.ndarray:
        """
        Create look-at view matrix.
        
        Args:
            eye: Camera position
            center: Point to look at
            up: Up vector
            
        Returns:
            4x4 view matrix
        """
        f = center - eye
        f = f / np.linalg.norm(f)
        
        s = np.cross(f, up)
        s = s / np.linalg.norm(s)
        
        u = np.cross(s, f)
        
        m = np.eye(4, dtype=np.float32)
        m[0, :3] = s
        m[1, :3] = u
        m[2, :3] = -f
        m[:3, 3] = -np.array([s, u, -f]) @ eye
        
        return m
    
    @staticmethod
    def perspective_projection(fov: float, aspect: float, near: float, far: float) -> np.ndarray:
        """
        Create perspective projection matrix.
        
        Args:
            fov: Field of view in radians
            aspect: Aspect ratio (width/height)
            near: Near plane distance
            far: Far plane distance
            
        Returns:
            4x4 projection matrix
        """
        f = 1.0 / math.tan(fov / 2.0)
        m = np.zeros((4, 4), dtype=np.float32)
        
        m[0, 0] = f / aspect
        m[1, 1] = f
        m[2, 2] = (far + near) / (near - far)
        m[2, 3] = (2.0 * far * near) / (near - far)
        m[3, 2] = -1.0
        
        return m


class OrbitCamera(Camera):
    """
    Orbiting camera that rotates around a target point.
    
    Suitable for general 3D exploration of storm structures.
    """
    
    def __init__(self, distance: float = 50.0, height: float = 20.0,
                 angle: float = 0.0, target: Optional[np.ndarray] = None):
        """
        Initialize orbit camera.
        
        Args:
            distance: Distance from target
            height: Height above ground
            angle: Azimuth angle in radians
            target: Target point to orbit around (default: origin)
        """
        self.distance = distance
        self.height = height
        self.angle = angle
        self.target = target if target is not None else np.array([0.0, 0.0, 10.0])
        self.up = np.array([0.0, 0.0, 1.0])
        
        # Projection parameters
        self.fov = math.radians(45.0)
        self.near = 0.1
        self.far = 200.0
    
    def set_angle(self, angle: float):
        """Set azimuth angle"""
        self.angle = angle
    
    def set_distance(self, distance: float):
        """Set distance from target"""
        self.distance = distance
    
    def set_height(self, height: float):
        """Set height above ground"""
        self.height = height
    
    def get_camera_position(self) -> np.ndarray:
        """Get camera position"""
        eye_x = self.distance * math.cos(self.angle)
        eye_y = self.distance * math.sin(self.angle)
        eye_z = self.height
        return np.array([eye_x, eye_y, eye_z])
    
    def get_view_matrix(self) -> np.ndarray:
        """Get view matrix"""
        eye = self.get_camera_position()
        return self.look_at(eye, self.target, self.up)
    
    def get_projection_matrix(self, width: int, height: int) -> np.ndarray:
        """Get projection matrix"""
        aspect = width / height if height > 0 else 1.0
        return self.perspective_projection(self.fov, aspect, self.near, self.far)


class ScientificCamera(Camera):
    """
    Scientific camera for side-view visualization of downdrafts.
    
    Provides a side-view angle suitable for visualizing descending columns
    and their interaction with the surface, similar to Lewellen et al. (2008).
    """
    
    def __init__(self, distance: float = 80.0, height: float = 15.0,
                 azimuth: float = math.pi / 2, elevation: float = 0.1,
                 target: Optional[np.ndarray] = None):
        """
        Initialize scientific camera.
        
        Args:
            distance: Distance from target
            height: Height offset
            azimuth: Azimuth angle (0 = front, π/2 = side)
            elevation: Elevation angle (0 = horizontal, positive = looking down)
            target: Target point (default: origin)
        """
        self.distance = distance
        self.height = height
        self.azimuth = azimuth
        self.elevation = elevation
        self.target = target if target is not None else np.array([0.0, 0.0, 5.0])
        self.up = np.array([0.0, 0.0, 1.0])
        
        # Projection parameters
        self.fov = math.radians(30.0)  # Narrower FOV for scientific view
        self.near = 0.1
        self.far = 300.0
    
    def get_camera_position(self) -> np.ndarray:
        """Get camera position"""
        # Calculate position based on azimuth and elevation
        eye_x = self.distance * math.cos(self.azimuth) * math.cos(self.elevation)
        eye_y = self.distance * math.sin(self.azimuth) * math.cos(self.elevation)
        eye_z = self.height + self.distance * math.sin(self.elevation)
        return np.array([eye_x, eye_y, eye_z])
    
    def get_view_matrix(self) -> np.ndarray:
        """Get view matrix"""
        eye = self.get_camera_position()
        return self.look_at(eye, self.target, self.up)
    
    def get_projection_matrix(self, width: int, height: int) -> np.ndarray:
        """Get projection matrix"""
        aspect = width / height if height > 0 else 1.0
        return self.perspective_projection(self.fov, aspect, self.near, self.far)
    
    def set_side_view(self):
        """Set to side view (azimuth = π/2)"""
        self.azimuth = math.pi / 2
        self.elevation = 0.1
    
    def set_front_view(self):
        """Set to front view (azimuth = 0)"""
        self.azimuth = 0.0
        self.elevation = 0.1
