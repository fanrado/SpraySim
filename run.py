#!/usr/bin/env python3
"""Command-line entry point for SpraySim.

Examples
--------
    python run.py                          # default spray, writes output/
    python run.py --droplets 8000 --speed 12 --cone 35
    python run.py --height 2.0 --out output/tall_spray.png
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from spraysim import SimConfig, NozzleConfig, PhysicsConfig, Simulator, analysis, plots


def build_config(args: argparse.Namespace) -> SimConfig:
    nozzle = NozzleConfig(
        position=(0.0, 0.0, args.height),
        half_angle=math.radians(args.cone),
        exit_speed=args.speed,
        mean_radius=args.radius_mm * 1e-3,
    )
    return SimConfig(
        n_droplets=args.droplets,
        dt=args.dt,
        seed=args.seed,
        nozzle=nozzle,
        physics=PhysicsConfig(),
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Physics-based spray droplet simulation.")
    p.add_argument("--droplets", type=int, default=4000, help="number of droplets")
    p.add_argument("--speed", type=float, default=9.0, help="nozzle exit speed (m/s)")
    p.add_argument("--cone", type=float, default=25.0, help="cone half-angle (degrees)")
    p.add_argument("--height", type=float, default=1.5, help="nozzle height (m)")
    p.add_argument("--radius-mm", type=float, default=0.4, help="mean droplet radius (mm)")
    p.add_argument("--dt", type=float, default=1e-3, help="integration timestep (s)")
    p.add_argument("--seed", type=int, default=42, help="RNG seed")
    p.add_argument("--out", type=Path, default=Path("output/spray_summary.png"),
                   help="output figure path")
    p.add_argument("--no-plot", action="store_true", help="skip figure generation")
    args = p.parse_args()

    config = build_config(args)

    print(f"Simulating {config.n_droplets} droplets "
          f"(speed={args.speed} m/s, cone={args.cone} deg, height={args.height} m)...")
    result = Simulator(config).run()
    stats = analysis.summarize(result, config)

    print("\n=== Spray statistics ===")
    print(json.dumps(stats.as_dict(), indent=2))

    if not args.no_plot:
        path = plots.save_figure(result, config, args.out)
        print(f"\nFigure written to {path.resolve()}")


if __name__ == "__main__":
    main()
