# Analysis

Offline analysis of saved SpraySim runs. Every simulation writes a
self-describing `.npz` archive to `output/` (per-droplet arrays **plus** the
config that produced them — see `spraysim/storage.py`), so these scripts can
regenerate every plot and statistic without re-running the simulation.

## `report.py` — PDF report

Builds a multi-page **PDF** from one or more `output/*.npz` archives.

```bash
# All runs in output/ -> output/spray_report.pdf
python analysis/report.py

# A single run
python analysis/report.py output/big_drops.npz

# Specific runs to a custom path
python analysis/report.py output/fine_mist.npz output/big_drops.npz --out report.pdf

# Custom glob
python analysis/report.py --glob "output/*.npz" --out output/spray_report.pdf
```

### What's in the PDF

- **Cover page** — when it was generated and the runs included.
- **Comparison table** *(only with >1 run)* — droplet count, exit speed, flow
  rate, mean flight time, p50/p90 coverage radius and mean radius side by side.
- **Per run** (three pages each):
  1. the standard 2×2 summary figure (trajectories, landing pattern, radial
     profile, size histogram);
  2. an **extra-analysis** page — radial coverage CDF and a droplet-size-vs-range
     scatter (does size sorting push bigger drops further?);
  3. a **configuration & statistics** page listing the exact material, nozzle,
     hydraulics and derived stats behind the run.

Archives that can't be read (e.g. an older, incompatible format, or a foreign
`.npz`) are skipped with a warning rather than aborting the whole report.

## `validate.py` — physics validation

Runs the simulator against closed-form benchmarks and prints a pass/fail table,
exiting non-zero if any check fails:

```bash
python analysis/validate.py
```

Checks: vacuum free-fall, terminal velocity (both the constant and Clift-Gauvin
drag models), the `Cd(Re)` Stokes/Newton limits, O(dt) timestep convergence, drag
monotonicity, the Torricelli / density-scaling hydraulics identities, impact-speed
energy consistency at the ground crossing, and the clipped-normal `E[r^3]`
correction. The cheap deterministic checks are also asserted in
`tests/test_simulator.py`; this script keeps the full battery (including the
slower terminal-velocity runs) in one runnable place. Takes a few seconds.

## Prerequisites

Uses `numpy` and `matplotlib` from `requirements.txt` — no extra dependencies.
Produce some input first by running a simulation (e.g. `./main.sh`), which writes
`output/*.npz`.
