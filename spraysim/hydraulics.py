"""Nozzle hydraulics: turn pressure + nozzle geometry into flow rate, exit
speed, and a droplet count.

Model (incompressible orifice flow, Bernoulli + Torricelli):

    ideal velocity   v_ideal = sqrt(2 * dP / rho)
    exit speed       v       = C_v * v_ideal
    volumetric flow  Q       = C_d * A * v_ideal          (A = orifice area)
    droplet count    N       = Q * spray_time / V_droplet

where C_d (discharge) and C_v (velocity) coefficients depend on nozzle shape,
and V_droplet is the mean droplet volume of the size distribution.
"""

from __future__ import annotations

import math

# Per-shape (discharge C_d, velocity C_v) coefficients. C_d folds in the
# vena-contracta contraction; C_v is the velocity loss. Representative values.
SHAPE_COEFFICIENTS: dict[str, tuple[float, float]] = {
    "sharp_orifice": (0.61, 0.98),    # thin, sharp-edged hole
    "rounded_orifice": (0.92, 0.98),  # well-rounded inlet, little contraction
    "full_cone": (0.75, 0.95),        # solid-cone agricultural nozzle
    "hollow_cone": (0.70, 0.95),      # swirl / hollow-cone atomiser
    "flat_fan": (0.88, 0.97),         # flat-fan spray tip
}

DEFAULT_SHAPE = "full_cone"


def shape_coefficients(shape: str) -> tuple[float, float]:
    """Return ``(C_d, C_v)`` for a named nozzle shape."""
    try:
        return SHAPE_COEFFICIENTS[shape]
    except KeyError:
        options = ", ".join(sorted(SHAPE_COEFFICIENTS))
        raise ValueError(
            f"Unknown nozzle shape {shape!r}. Choose one of: {options}."
        ) from None


def orifice_area(diameter: float) -> float:
    """Cross-sectional area (m^2) of a circular orifice of given diameter (m)."""
    return math.pi * (diameter / 2.0) ** 2


def ideal_velocity(pressure: float, density: float) -> float:
    """Torricelli ideal jet velocity (m/s) for a pressure drop ``pressure`` (Pa)."""
    if pressure <= 0.0:
        return 0.0
    return math.sqrt(2.0 * pressure / density)


def exit_speed(pressure: float, shape: str, density: float) -> float:
    """Actual droplet launch speed (m/s): ``C_v`` times the ideal velocity."""
    _, cv = shape_coefficients(shape)
    return cv * ideal_velocity(pressure, density)


def flow_rate(pressure: float, diameter: float, shape: str, density: float) -> float:
    """Volumetric flow rate (m^3/s) through the orifice."""
    cd, _ = shape_coefficients(shape)
    return cd * orifice_area(diameter) * ideal_velocity(pressure, density)


def droplet_count(
    flow: float,
    spray_duration: float,
    mean_droplet_volume: float,
    max_droplets: int | None = None,
) -> int:
    """Number of droplets produced: total sprayed volume / mean droplet volume.

    Always at least 1 (so a valid run is possible), and capped at
    ``max_droplets`` when given.
    """
    if mean_droplet_volume <= 0.0:
        raise ValueError("mean_droplet_volume must be positive.")
    n = int(round(flow * spray_duration / mean_droplet_volume))
    n = max(1, n)
    if max_droplets is not None:
        n = min(n, max_droplets)
    return n
