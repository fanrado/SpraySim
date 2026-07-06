#!/usr/bin/env python3
"""Physics validation harness for SpraySim.

Runs the simulator against closed-form benchmarks and prints a pass/fail table.
Exits non-zero if any check fails, so it doubles as a smoke test for the physics.
The cheap, deterministic checks are also asserted in ``tests/test_simulator.py``;
this script keeps the full battery (including the slower terminal-velocity runs)
in one runnable place.

    python analysis/validate.py            # run all checks
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

# Allow "python analysis/validate.py" from the repo root without installing.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from spraysim import (  # noqa: E402
    SimConfig, NozzleConfig, PhysicsConfig, PathConfig, Simulator,
    analysis, drag, hydraulics,
)
from spraysim.nozzle import Nozzle  # noqa: E402

G = 9.81
RHO_AIR = 1.225
MU_AIR = drag.DEFAULT_AIR_VISCOSITY
WATER = 1000.0

_results: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    _results.append((name, bool(ok), detail))


def _drop_from_rest(r, model, height, dt, air=RHO_AIR, max_time=60.0):
    """A single droplet of fixed radius released from rest, integrated to impact."""
    noz = NozzleConfig(position=(0, 0, height), direction=(0, 0, -1), half_angle=0.0,
                       pressure=0.0, speed_spread=0.0,
                       distribution="normal", mean_radius=r, radius_std=0.0)
    cfg = SimConfig(n_droplets=1, dt=dt, max_time=max_time, seed=0, nozzle=noz,
                    physics=PhysicsConfig(air_density=air, drag_model=model))
    return Simulator(cfg).run()


def _terminal_cg(r):
    """Fixed-point solve of g = k(Cd(Re)) v^2 with the Clift-Gauvin Cd(Re)."""
    v = 5.0
    for _ in range(300):
        Re = RHO_AIR * v * 2 * r / MU_AIR
        cd = float(drag.clift_gauvin(np.array(Re)))
        v = math.sqrt(G / (3.0 * RHO_AIR * cd / (8.0 * WATER * r)))
    return v


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #

def check_vacuum_freefall():
    h = 20.0
    res = _drop_from_rest(4e-4, "clift_gauvin", h, 1e-4, air=0.0, max_time=10.0)
    t_an, v_an = math.sqrt(2 * h / G), G * math.sqrt(2 * h / G)
    et = abs(res.flight_times[0] - t_an) / t_an
    ev = abs(res.impact_speeds[0] - v_an) / v_an
    record("vacuum free-fall: t=sqrt(2h/g), v=g t", et < 1e-3 and ev < 1e-3,
           f"t err {et:.1e}, v err {ev:.1e}")


def check_terminal_velocity_constant():
    ok, parts = True, []
    for r in (1e-4, 4e-4, 1e-3):
        vt = math.sqrt(G / (3.0 * RHO_AIR * 0.47 / (8.0 * WATER * r)))
        res = _drop_from_rest(r, "constant", height=40.0, dt=1e-4)
        e = abs(res.impact_speeds[0] - vt) / vt
        ok &= e < 5e-3
        parts.append(f"{r*1e3:.1f}mm:{e:.1e}")
    record("terminal velocity = sqrt(g/k) [constant Cd]", ok, ", ".join(parts))


def check_terminal_velocity_reynolds():
    ok, parts = True, []
    for r in (4e-4, 8e-4, 1.6e-3):
        vt = _terminal_cg(r)
        res = _drop_from_rest(r, "clift_gauvin", height=40.0, dt=1e-4)
        e = abs(res.impact_speeds[0] - vt) / vt
        ok &= e < 1e-2
        parts.append(f"{r*1e3:.1f}mm:{e:.1e}")
    record("terminal velocity = Cd(Re) fixed point [clift_gauvin]", ok, ", ".join(parts))


def check_drag_coefficient_limits():
    lo = float(drag.clift_gauvin(np.array(0.01)))
    hi = float(drag.clift_gauvin(np.array(1e5)))
    ok = abs(lo - 24.0 / 0.01) / (24.0 / 0.01) < 0.02 and 0.35 < hi < 0.55
    record("Cd(Re): Stokes 24/Re at low Re, Newton plateau at high Re", ok,
           f"Cd(0.01)={lo:.1f}~{24/0.01:.0f}, Cd(1e5)={hi:.3f}")


def check_timestep_convergence():
    h = 10.0
    ref = _drop_from_rest(4e-4, "constant", h, 1.25e-4, max_time=30.0).flight_times[0]
    errs = [abs(_drop_from_rest(4e-4, "constant", h, dt, max_time=30.0).flight_times[0] - ref)
            for dt in (4e-3, 2e-3, 1e-3)]
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    ok = all(1.6 < rt < 2.4 for rt in ratios)
    record("timestep convergence is O(dt) (error halves as dt halves)", ok,
           "ratios " + ", ".join(f"{rt:.2f}" for rt in ratios))


def check_drag_monotonicity():
    ranges = []
    for r in (5e-5, 1e-4, 2e-4, 4e-4, 8e-4, 1.6e-3):
        noz = NozzleConfig(position=(0, 0, 2.0), direction=(1, 0, 0), half_angle=0.0,
                           pressure=3e5, speed_spread=0.0,
                           distribution="normal", mean_radius=r, radius_std=0.0)
        res = Simulator(SimConfig(n_droplets=1, dt=1e-3, max_time=20.0, seed=0,
                                  nozzle=noz)).run()
        ranges.append(math.hypot(res.landing_positions[0, 0], res.landing_positions[0, 1]))
    ok = all(ranges[i] < ranges[i + 1] for i in range(len(ranges) - 1))
    record("smaller droplets travel less far (drag monotonicity)", ok,
           "ranges(m) " + ", ".join(f"{x:.2f}" for x in ranges))


def check_hydraulics_identities():
    ok = all(abs(hydraulics.ideal_velocity(P, WATER) - math.sqrt(2 * P / WATER)) < 1e-9
             for P in (1e5, 3e5, 5e5))
    vw = hydraulics.exit_speed(3e5, "full_cone", 1000.0)
    ve = hydraulics.exit_speed(3e5, "full_cone", 789.0)
    ok &= abs(ve / vw - math.sqrt(1000.0 / 789.0)) < 1e-6
    record("hydraulics: Torricelli & 1/sqrt(rho) density scaling", ok,
           f"v(789)/v(1000)={ve/vw:.5f}")


def check_impact_speed_energy():
    h = 2.0
    noz = NozzleConfig(position=(0, 0, h), direction=(1, 0, 0), half_angle=0.0,
                       pressure=5e5, speed_spread=0.0,
                       distribution="normal", mean_radius=8e-4, radius_std=0.0)
    res = Simulator(SimConfig(n_droplets=1, dt=1e-3, max_time=20.0, seed=0, nozzle=noz,
                              physics=PhysicsConfig(air_density=0.0))).run()
    u, v = res.launch_speeds[0], res.impact_speeds[0]
    v_energy = math.sqrt(u * u + 2 * G * h)
    e = abs(v - v_energy) / v_energy
    record("impact speed energy-consistent at crossing (vacuum)", e < 5e-4,
           f"err {e:.1e}")


def check_normal_moment():
    rng = np.random.default_rng(0)
    ok, parts = True, []
    for m, s in ((4e-4, 1e-4), (2e-4, 3e-4)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nz = Nozzle(NozzleConfig(distribution="normal", mean_radius=m, radius_std=s), WATER)
            samp = float(np.mean(nz.sample_radii(500_000, rng) ** 3))
        e = abs(nz.mean_cubed_radius() - samp) / samp
        ok &= e < 1e-2
        parts.append(f"s/m={s/m:.1f}:{e:.1e}")
    record("normal E[r^3] matches clipped sampler (any spread)", ok, ", ".join(parts))


def _raster_gcode(size_mm=200, pitch_mm=10, standoff_mm=150, feed=3000):
    lines = ["G21", "G90", f"G0 X0 Y0 Z{standoff_mm}", f"G1 F{feed}"]
    rows = size_mm // pitch_mm + 1
    for i in range(rows):
        y = i * pitch_mm
        x = size_mm if i % 2 == 0 else 0
        lines.append(f"G1 X{x} Y{y}")
        if i < rows - 1:
            lines.append(f"G0 Y{y + pitch_mm}")
    return "\n".join(lines) + "\n"


def check_path_uniformity():
    """A raster toolpath builds a far more uniform interior film than one spot."""
    roi = (0.04, 0.16, 0.04, 0.16)
    cfg = SimConfig(n_droplets=15000, seed=3, path=PathConfig(gcode=_raster_gcode()))
    raster = analysis.uniformity(
        analysis.deposition_map(Simulator(cfg).run(), cfg, cell_size=0.02), roi=roi)
    spot = analysis.uniformity(
        analysis.deposition_map(Simulator(SimConfig(n_droplets=15000, seed=3)).run(),
                                SimConfig(), cell_size=0.02), roi=roi)
    ok = raster.christiansen_cu > 0.7 and raster.christiansen_cu > spot.christiansen_cu
    record("raster path builds a uniform film (CU) vs a single spot", ok,
           f"raster CU={raster.christiansen_cu:.2f} vs spot CU={spot.christiansen_cu:.2f}")


CHECKS = [
    check_vacuum_freefall,
    check_terminal_velocity_constant,
    check_terminal_velocity_reynolds,
    check_drag_coefficient_limits,
    check_timestep_convergence,
    check_drag_monotonicity,
    check_hydraulics_identities,
    check_impact_speed_energy,
    check_normal_moment,
    check_path_uniformity,
]


def main() -> int:
    print("SpraySim physics validation\n" + "=" * 60)
    for fn in CHECKS:
        fn()
    for name, ok, detail in _results:
        tag = "PASS" if ok else "FAIL"
        print(f"[{tag}] {name}")
        if detail:
            print(f"       {detail}")
    n_fail = sum(not ok for _, ok, _ in _results)
    print("=" * 60)
    print(f"{len(_results) - n_fail}/{len(_results)} checks passed"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
