"""Droplet drag coefficient models.

The aerodynamic drag on a droplet is ``F = 0.5 * rho_air * Cd * A * |v| * v``.
The drag coefficient ``Cd`` is not really constant: it depends on the droplet
Reynolds number ``Re = rho_air * |v| * d / mu_air`` (d = 2r). This module supplies
``Cd(Re)`` for the available models.

Models
------
``clift_gauvin`` (default)
    Smooth correlation valid across the full sub-critical range (Re up to ~2e5):

        Cd = (24/Re)*(1 + 0.15*Re^0.687) + 0.42 / (1 + 4.25e4*Re^-1.16)

    Reduces to Stokes drag (Cd = 24/Re) as Re -> 0 and to ~0.44 (Newton) at high
    Re, so small droplets are slowed correctly while large droplets match the old
    behaviour.

``constant``
    A fixed ``Cd`` (``PhysicsConfig.drag_coefficient``, default 0.47). This is the
    original model; kept for reproducing pre-existing runs and for A/B tests.
"""

from __future__ import annotations

import numpy as np

CLIFT_GAUVIN = "clift_gauvin"
CONSTANT = "constant"

DRAG_MODELS = (CLIFT_GAUVIN, CONSTANT)

DEFAULT_DRAG_MODEL = CLIFT_GAUVIN
DEFAULT_AIR_VISCOSITY = 1.81e-5  # Pa*s, dry air at ~15 C

# Below this Reynolds number the drag force is negligible anyway; clamp to keep
# the 24/Re term finite (at v=0 the |v|*v factor is zero, so Cd*0 -> 0 cleanly).
_RE_FLOOR = 1e-12


def reynolds_number(speed, radius, air_density: float, air_viscosity: float):
    """Droplet Reynolds number ``Re = rho_air * speed * (2 r) / mu_air``."""
    return air_density * speed * (2.0 * radius) / air_viscosity


def clift_gauvin(Re):
    """Clift-Gauvin drag coefficient (array-safe)."""
    Re = np.maximum(Re, _RE_FLOOR)
    return (24.0 / Re) * (1.0 + 0.15 * Re**0.687) + 0.42 / (1.0 + 4.25e4 * Re**-1.16)


def drag_coefficient(Re, model: str = DEFAULT_DRAG_MODEL, constant: float = 0.47):
    """Return ``Cd`` for the given ``Re`` and model.

    ``Re`` may be a scalar or an array. For ``constant`` the ``Re`` argument is
    ignored and ``constant`` is returned (broadcast to ``Re``'s shape).
    """
    if model == CLIFT_GAUVIN:
        return clift_gauvin(Re)
    if model == CONSTANT:
        return np.broadcast_to(np.asarray(constant, dtype=float), np.shape(Re))
    raise ValueError(
        f"Unknown drag model {model!r}; choose one of: {', '.join(DRAG_MODELS)}."
    )
