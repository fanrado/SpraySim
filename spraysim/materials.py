"""Sprayed-liquid material properties.

The simulation is no longer water-specific: the sprayed liquid's **density**
(kg/m^3) feeds both the droplet mass in the drag term and the Torricelli exit
speed / flow rate in :mod:`spraysim.hydraulics`. This module holds a small
registry of common liquids so a material can be picked by name, with the density
overridable for anything not listed.

Densities are nominal values near 15-20 C.
"""

from __future__ import annotations

# name -> density (kg/m^3)
MATERIALS: dict[str, float] = {
    "water": 1000.0,
    "seawater": 1025.0,
    "ethanol": 789.0,
    "methanol": 792.0,
    "acetone": 784.0,
    "gasoline": 745.0,
    "kerosene": 810.0,
    "diesel": 832.0,
    "olive_oil": 915.0,
    "glycerin": 1260.0,
}

DEFAULT_MATERIAL = "water"


def material_density(name: str) -> float:
    """Density (kg/m^3) of a named material.

    Raises ``KeyError`` with the list of known materials if ``name`` is unknown.
    """
    try:
        return MATERIALS[name]
    except KeyError:
        known = ", ".join(sorted(MATERIALS))
        raise KeyError(f"Unknown material {name!r}; known materials: {known}") from None
