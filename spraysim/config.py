"""Configuration dataclasses for a spray simulation.

All quantities are SI units (metres, seconds, kilograms, radians).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math


@dataclass
class PhysicsConfig:
    """Environmental / fluid constants used by the integrator."""

    gravity: float = 9.81           # m/s^2, acts along -z
    air_density: float = 1.225      # kg/m^3 at sea level, 15 C
    water_density: float = 1000.0   # kg/m^3, droplet material
    drag_coefficient: float = 0.47  # dimensionless, sphere
    ground_z: float = 0.0           # m, height of the impact plane


@dataclass
class NozzleConfig:
    """Describes where droplets are born and how they are launched."""

    # Nozzle location (m). Default: 1.5 m above the origin.
    position: tuple[float, float, float] = (0.0, 0.0, 1.5)
    # Spray axis direction (need not be normalised). Default: straight down.
    direction: tuple[float, float, float] = (0.0, 0.0, -1.0)
    # Full cone spread is 2 * half_angle around the axis.
    half_angle: float = math.radians(25.0)
    # Exit speed (m/s) and its relative standard deviation.
    exit_speed: float = 9.0
    speed_spread: float = 0.15
    # Droplet radius distribution (m), sampled log-normally.
    mean_radius: float = 4.0e-4      # 0.4 mm
    radius_spread: float = 0.35      # relative sigma of the log-normal


@dataclass
class SimConfig:
    """Top-level run configuration."""

    n_droplets: int = 4000
    dt: float = 1.0e-3               # s, integration timestep
    max_time: float = 8.0            # s, hard cap on flight time
    n_trajectories: int = 60         # droplets whose full path is recorded
    seed: int | None = 42

    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    nozzle: NozzleConfig = field(default_factory=NozzleConfig)
