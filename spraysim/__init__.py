"""SpraySim — a physics-based spray/droplet particle simulator.

A nozzle emits droplets into a velocity cone. Each droplet is integrated
under gravity and aerodynamic drag until it lands on the ground plane.
The package produces numerical results and static plots (no animation).
"""

from .config import PhysicsConfig, NozzleConfig, SimConfig
from .nozzle import Nozzle
from .simulator import Simulator, SimResult
from . import analysis, plots

__all__ = [
    "PhysicsConfig",
    "NozzleConfig",
    "SimConfig",
    "Nozzle",
    "Simulator",
    "SimResult",
    "analysis",
    "plots",
]

__version__ = "0.1.0"
