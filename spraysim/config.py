"""Configuration dataclasses for a spray simulation.

All quantities are SI units (metres, seconds, kilograms, pascals, radians).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

from .hydraulics import DEFAULT_SHAPE
from .materials import DEFAULT_MATERIAL, DEFAULT_VISCOSITY, material_density
from .drag import DEFAULT_DRAG_MODEL, DEFAULT_AIR_VISCOSITY


@dataclass
class PhysicsConfig:
    """Environmental constants used by the integrator (the surrounding air)."""

    gravity: float = 9.81           # m/s^2, acts along -z
    air_density: float = 1.225      # kg/m^3 at sea level, 15 C
    air_viscosity: float = DEFAULT_AIR_VISCOSITY  # Pa*s; sets the droplet Reynolds number
    drag_coefficient: float = 0.47  # dimensionless; used only by the "constant" drag model
    drag_model: str = DEFAULT_DRAG_MODEL  # "clift_gauvin" (Re-dependent) | "constant"
    ground_z: float = 0.0           # m, height of the impact plane


@dataclass
class MaterialConfig:
    """The sprayed liquid. Its density sets droplet mass and exit hydraulics.

    ``viscosity`` (dynamic, Pa*s) is carried for custom-liquid definitions and
    reporting; it defaults to water's value and is not yet coupled to the flight
    physics.

    ``solids_fraction`` is the volume fraction of solid material in the prepared
    solution (the rest is solvent that evaporates). It sets the *dry* deposited
    film thickness = wet thickness * solids_fraction; ``1.0`` = pure liquid.
    """

    name: str = DEFAULT_MATERIAL              # label; picks default properties
    density: float = material_density(DEFAULT_MATERIAL)  # kg/m^3
    viscosity: float = DEFAULT_VISCOSITY      # Pa*s (dynamic); water by default
    solids_fraction: float = 1.0              # volume fraction of solids (dry film)


@dataclass
class NozzleConfig:
    """Describes where droplets are born and how they are launched.

    Exit speed, flow rate and droplet count are *not* set here — they are
    derived from ``pressure``, ``orifice_diameter`` and ``shape`` by the
    hydraulics model (see :mod:`spraysim.hydraulics`).
    """

    # Nozzle location (m). Default: 1.5 m above the origin.
    position: tuple[float, float, float] = (0.0, 0.0, 1.5)
    # Spray axis direction (need not be normalised). Default: straight down.
    direction: tuple[float, float, float] = (0.0, 0.0, -1.0)
    # Full cone spread is 2 * half_angle around the axis.
    half_angle: float = math.radians(25.0)

    # --- Hydraulics inputs (these drive exit speed, flow rate, droplet count) ---
    pressure: float = 3.0e5             # Pa (3 bar) gauge pressure drop
    orifice_diameter: float = 0.8e-3    # m (0.8 mm) orifice diameter
    shape: str = DEFAULT_SHAPE          # nozzle shape -> discharge/velocity coeffs
    # Relative std of droplet speed about the derived exit speed (turbulence).
    speed_spread: float = 0.15

    # --- Droplet size distribution ---
    # "normal" (Gaussian) or "lognormal"; both parameterised by the linear-space
    # mean radius and standard deviation below.
    distribution: str = "lognormal"
    mean_radius: float = 4.0e-4         # m (0.4 mm) mean droplet radius
    radius_std: float = 1.2e-4          # m (0.12 mm) std of droplet radius


@dataclass
class SimConfig:
    """Top-level run configuration."""

    # Droplets are normally derived from hydraulics. Set this to pin an explicit
    # count (used by tests and quick experiments) and skip the derivation.
    n_droplets: int | None = None
    spray_duration: float = 0.15        # s the nozzle is open (sets droplet count)
    max_droplets: int = 200_000         # safety cap on the derived droplet count

    dt: float = 1.0e-3                  # s, integration timestep
    max_time: float = 8.0               # s, hard cap on flight time
    n_trajectories: int = 60            # droplets whose full path is recorded
    seed: int | None = 42

    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    material: MaterialConfig = field(default_factory=MaterialConfig)
    nozzle: NozzleConfig = field(default_factory=NozzleConfig)
