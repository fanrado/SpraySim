#!/usr/bin/env python3
"""Command-line entry point for SpraySim.

Droplet count and exit speed are derived from the nozzle pressure, orifice size
and shape; the droplet size distribution is configurable (normal / lognormal).

Examples
--------
    python run.py                                   # default nozzle, writes output/
    python run.py --pressure-bar 5 --orifice-mm 1.2 --shape flat_fan
    python run.py --distribution normal --mean-radius-mm 0.3 --radius-std-mm 0.08
    python run.py --droplets 5000                   # pin an explicit droplet count
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from spraysim import (
    SimConfig,
    NozzleConfig,
    PhysicsConfig,
    MaterialConfig,
    Simulator,
    analysis,
    plots,
    storage,
)
from spraysim.hydraulics import SHAPE_COEFFICIENTS
from spraysim.materials import (
    MATERIALS,
    DEFAULT_VISCOSITY,
    material_density,
    material_viscosity,
)
from spraysim.drag import DRAG_MODELS, DEFAULT_DRAG_MODEL


def build_config(args: argparse.Namespace) -> SimConfig:
    # Material density: explicit --density overrides the named material's default.
    density = args.density if args.density is not None else material_density(args.material)
    # Material viscosity: --viscosity wins; else the registry value; else (custom
    # liquid not in the registry) fall back to water's viscosity.
    if args.viscosity is not None:
        viscosity = args.viscosity
    elif args.material in MATERIALS:
        viscosity = material_viscosity(args.material)
    else:
        viscosity = DEFAULT_VISCOSITY
    material = MaterialConfig(name=args.material, density=density, viscosity=viscosity,
                              solids_fraction=args.solids_fraction)
    nozzle = NozzleConfig(
        position=(0.0, 0.0, args.height),
        half_angle=math.radians(args.cone),
        pressure=args.pressure_bar * 1e5,          # bar -> Pa
        orifice_diameter=args.orifice_mm * 1e-3,   # mm -> m
        shape=args.shape,
        speed_spread=args.speed_spread,
        distribution=args.distribution,
        mean_radius=args.mean_radius_mm * 1e-3,    # mm -> m
        radius_std=args.radius_std_mm * 1e-3,      # mm -> m
    )
    return SimConfig(
        n_droplets=args.droplets,                  # None => derive from hydraulics
        spray_duration=args.spray_duration,
        dt=args.dt,
        seed=args.seed,
        nozzle=nozzle,
        material=material,
        physics=PhysicsConfig(drag_model=args.drag_model),
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Physics-based spray droplet simulation.")
    # Hydraulics (drive exit speed, flow rate and droplet count).
    p.add_argument("--pressure-bar", type=float, default=3.0, help="nozzle pressure (bar)")
    p.add_argument("--orifice-mm", type=float, default=0.8, help="orifice diameter (mm)")
    p.add_argument("--shape", default="full_cone", choices=sorted(SHAPE_COEFFICIENTS),
                   help="nozzle shape (sets discharge/velocity coefficients)")
    p.add_argument("--spray-duration", type=float, default=0.15,
                   help="how long the nozzle is open (s); scales droplet count")
    p.add_argument("--droplets", type=int, default=None,
                   help="override the derived droplet count with an explicit number")
    # Sprayed liquid material (density drives droplet mass and exit hydraulics).
    p.add_argument("--material", default="water",
                   help="sprayed liquid; known: " + ", ".join(sorted(MATERIALS))
                        + " (or any name with --density)")
    p.add_argument("--density", type=float, default=None,
                   help="liquid density (kg/m^3); overrides the material's default")
    p.add_argument("--viscosity", type=float, default=None,
                   help="liquid dynamic viscosity (Pa*s); overrides the material's "
                        "default (custom liquids default to water's viscosity)")
    p.add_argument("--solids-fraction", type=float, default=1.0,
                   help="volume fraction of solids in the sprayed solution "
                        "(sets dry film thickness; 1.0 = pure liquid)")
    # Droplet size distribution.
    p.add_argument("--distribution", default="lognormal", choices=["normal", "lognormal"],
                   help="droplet radius distribution")
    p.add_argument("--mean-radius-mm", type=float, default=0.4, help="mean droplet radius (mm)")
    p.add_argument("--radius-std-mm", type=float, default=0.12,
                   help="std of droplet radius (mm)")
    # Physics.
    p.add_argument("--drag-model", default=DEFAULT_DRAG_MODEL, choices=list(DRAG_MODELS),
                   help="droplet drag model: clift_gauvin (Reynolds-dependent) or "
                        "constant (fixed Cd; legacy behaviour)")
    # Geometry / integration.
    p.add_argument("--cone", type=float, default=25.0, help="cone half-angle (degrees)")
    p.add_argument("--height", type=float, default=1.5, help="nozzle height (m)")
    p.add_argument("--speed-spread", type=float, default=0.15,
                   help="relative std of droplet speed about the exit speed")
    p.add_argument("--dt", type=float, default=1e-3, help="integration timestep (s)")
    p.add_argument("--seed", type=int, default=42, help="RNG seed")
    p.add_argument("--out", type=Path, default=Path("output/spray_summary.png"),
                   help="output figure path")
    p.add_argument("--no-plot", action="store_true", help="skip figure generation")
    p.add_argument("--data", type=Path, default=Path("output/spray_data.npz"),
                   help="output .npz data path (result arrays + config)")
    p.add_argument("--no-data", action="store_true", help="skip saving the .npz data")
    args = p.parse_args()

    try:
        config = build_config(args)
    except KeyError as exc:
        p.error(str(exc).strip('"'))

    print(f"Material: {config.material.name} ({config.material.density:g} kg/m^3, "
          f"{config.material.viscosity:g} Pa*s, solids {config.material.solids_fraction:g})")
    print(f"Nozzle: {args.pressure_bar} bar, orifice {args.orifice_mm} mm, "
          f"{args.shape}, spraying for {args.spray_duration} s")
    print(f"Droplet size: {args.distribution}, "
          f"mean {args.mean_radius_mm} mm +/- {args.radius_std_mm} mm")
    print(f"Drag model: {config.physics.drag_model}\n")

    result = Simulator(config).run()
    stats = analysis.summarize(result, config)

    print("=== Derived hydraulics ===")
    print(json.dumps({
        "exit_speed_ms": round(result.exit_speed, 3),
        "flow_rate_ml_s": round(result.flow_rate * 1e6, 3),
        "droplets": result.n,
        "droplets_capped": result.droplets_capped,
    }, indent=2))

    print("\n=== Spray statistics ===")
    print(json.dumps(stats.as_dict(), indent=2))

    if result.droplets_capped:
        print(f"\nNote: droplet count was capped at {config.max_droplets}. "
              "Lower --spray-duration or raise the cap for a full run.")

    if not args.no_data:
        data_path = storage.save_result(result, config, args.data)
        print(f"\nData written to {data_path.resolve()}")

    if not args.no_plot:
        path = plots.save_figure(result, config, args.out)
        print(f"Figure written to {path.resolve()}")


if __name__ == "__main__":
    main()
