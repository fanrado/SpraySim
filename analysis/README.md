# Analysis

Offline analysis of saved SpraySim runs. Every simulation writes a
self-describing `.npz` archive to its run folder `output/<run-name>/`
(`output/tmp/` for unnamed runs; per-droplet arrays **plus** the
config that produced them — see `spraysim/storage.py`), so these scripts can
regenerate every plot and statistic without re-running the simulation.

## Deposition & uniformity (library)

`spraysim.analysis` turns a run into a **dry film-thickness map** and
**coating-uniformity** metrics:

```python
from spraysim import storage, analysis
result, config = storage.load_result("output/tmp/spray_data.npz")

field = analysis.deposition_map(result, config)      # DepositionField (dry thickness, m)
stats = analysis.uniformity(field)                   # over the wetted area, or roi=(x0,x1,y0,y1)
print(stats.as_dict())
```

- `deposition_map(result, config, *, cell_size=None, extent=None)` bins landed
  droplet volumes onto a grid; the dry thickness folds in
  `material.solids_fraction`.
- `uniformity(field, *, roi=None, coverage_threshold=None)` returns the CV,
  Christiansen uniformity coefficient (CU), coverage fraction and thickness
  percentiles. Default ROI is the wetted region; pass a rectangle to score a
  target area (dry gaps penalised). A single-spot spray is highly non-uniform
  (CV ≈ 1, CU ≈ 0.3) — path spraying (G-code) is how you build a uniform film.

## `report.py` — PDF report

Builds a multi-page **PDF** from one or more `output/<run>/*.npz` archives.

```bash
# Every run folder under output/ -> output/spray_report.pdf
python analysis/report.py

# A single run
python analysis/report.py output/tmp/big_drops.npz

# Specific runs to a custom path
python analysis/report.py output/mist/fine_mist.npz output/drops/big_drops.npz --out report.pdf

# Custom glob (** recurses into run folders)
python analysis/report.py --glob "output/trial_*/*.npz" --out output/trials.pdf
```

### What's in the PDF

- **Cover page** — when it was generated and the runs included.
- **Comparison table** *(only with >1 run)* — droplet count, exit speed, flow
  rate, mean flight time, p90 coverage radius, mean radius and the Christiansen
  uniformity coefficient (CU) side by side.
- **Per run** (four pages each):
  1. the standard 2×2 summary figure (trajectories, landing pattern, radial
     profile, size histogram);
  2. an **extra-analysis** page — radial coverage CDF and a droplet-size-vs-range
     scatter (does size sorting push bigger drops further?);
  3. a **deposition & uniformity** page — dry film-thickness heatmap plus a
     per-cell thickness histogram annotated with CU / CV / coverage;
  4. a **configuration & statistics** page listing the exact material, nozzle,
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

## `fit_ptp_efficiency.py` — fit a WLS-efficiency measurement with a coating model

Fits the IFIC VUV measurement of p-terphenyl (pTP) sprayed from toluene onto
3 mm acrylic (re-emitted / incident photons vs. estimated pTP areal density) at
one wavelength, using the **simulated coating** as the model's morphology:

1. SpraySim runs `config/ptp_toluene.conf` (toluene + pTP, the real 120×120 mm
   raster in `examples/raster_120mm_pitch12.gcode`) and the dry-thickness map
   over the sample is normalised to mean 1 — the shape comes from the
   simulation, the amount of pTP from each sample's estimated density.
2. Locally the efficiency saturates, `ε(t) = ε_sat · (1 − exp(−t/t₀))`; the
   measurement is the area average over the sample. `t₀` (quoted as
   `ρ₀ = ρ_pTP · t₀` in µg/cm²) lumps the VUV absorption depth with the
   droplet-scale patchiness of the dried film (a Poisson overlap of dried
   droplet spots has exactly this form with a renormalised `ρ₀`).
3. `ε_sat` and `ρ₀` are fitted (weighted by the quoted SEM plus an assumed
   relative uncertainty on the estimated densities, `--rho-frac-err`), for each
   simulated spray-cone half-angle and for the uniform-film limit.

```bash
python analysis/fit_ptp_efficiency.py                        # 130 nm, six plotted samples
python analysis/fit_ptp_efficiency.py --cone-deg 5 10 15 25  # scan the cone angle
python analysis/fit_ptp_efficiency.py --samples all --wavelength 140
python analysis/fit_ptp_efficiency.py --run-name ptp_130nm   # -> output/ptp_130nm/
```

Writes `ptp_fit_130nm.png` into the run folder (`output/tmp/` unless
`--run-name` is given; data + fitted curves, the simulated thickness
map, its distribution) and a `.json` with every fit result and simulation
statistic. The data CSV path defaults to the IFIC QE-corrected export
(`--csv`). Needs `scipy` in addition to the base requirements.

## Prerequisites

Uses `numpy` and `matplotlib` from `requirements.txt` (`fit_ptp_efficiency.py`
also needs `scipy`).
Produce some input first by running a simulation (e.g. `./main.sh`), which writes
`output/<run-name>/*.npz`.
