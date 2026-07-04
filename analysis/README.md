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

## Prerequisites

Uses `numpy` and `matplotlib` from `requirements.txt` — no extra dependencies.
Produce some input first by running a simulation (e.g. `./main.sh`), which writes
`output/*.npz`.
