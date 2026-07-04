"""Sprayed-liquid material properties.

The simulation is not water-specific: the sprayed liquid is described by a couple
of physical properties held in a small registry so a material can be picked by
name, with values overridable for anything not listed.

* **density** (kg/m^3) feeds the droplet mass in the drag term and the Torricelli
  exit speed / flow rate in :mod:`spraysim.hydraulics`.
* **viscosity** (dynamic, Pa*s) is carried on the material for custom-liquid
  definitions and reporting. It is not yet coupled to the flight physics (which
  uses the surrounding *air*, not the liquid); it defaults to water's value.

Values are nominal at ~20 C.
"""

from __future__ import annotations

from typing import NamedTuple


class MaterialProps(NamedTuple):
    """Physical properties of a sprayed liquid."""

    density: float      # kg/m^3
    viscosity: float    # Pa*s (dynamic)


# name -> (density kg/m^3, dynamic viscosity Pa*s), nominal at ~20 C.
MATERIALS: dict[str, MaterialProps] = {
    "water":     MaterialProps(1000.0, 1.00e-3),
    "seawater":  MaterialProps(1025.0, 1.08e-3),
    "ethanol":   MaterialProps(789.0, 1.20e-3),
    "methanol":  MaterialProps(792.0, 0.59e-3),
    "acetone":   MaterialProps(784.0, 0.32e-3),
    "gasoline":  MaterialProps(745.0, 0.60e-3),
    "kerosene":  MaterialProps(810.0, 1.50e-3),
    "diesel":    MaterialProps(832.0, 2.50e-3),
    "olive_oil": MaterialProps(915.0, 84.0e-3),
    "glycerin":  MaterialProps(1260.0, 1.41),
}

DEFAULT_MATERIAL = "water"


def _lookup(name: str) -> MaterialProps:
    try:
        return MATERIALS[name]
    except KeyError:
        known = ", ".join(sorted(MATERIALS))
        raise KeyError(f"Unknown material {name!r}; known materials: {known}") from None


def material_density(name: str) -> float:
    """Density (kg/m^3) of a named material (raises KeyError if unknown)."""
    return _lookup(name).density


def material_viscosity(name: str) -> float:
    """Dynamic viscosity (Pa*s) of a named material (raises KeyError if unknown)."""
    return _lookup(name).viscosity


# Water's viscosity, used as the default for custom / unspecified liquids.
DEFAULT_VISCOSITY = MATERIALS[DEFAULT_MATERIAL].viscosity
