# SpraySim — Feature Plan: Deposition Analysis & G-code Path Spraying

Two feature areas:

- **A — Deposition analysis:** measure the **thickness** of the deposited film and
  the **uniformity** of the coating.
- **B — Path spraying:** move the nozzle along a toolpath written in **G-code**
  instead of spraying from a single fixed position.

They are complementary: the thickness/uniformity metrics (A) are exactly what you
use to judge the quality of a painted path (B). Build **A first** — it also gives
B a quantitative acceptance test (a dense raster should yield a uniform film).

Current state this builds on:
- `Simulator.run()` integrates a *batch* of independent droplets and returns a
  `SimResult` with `landing_positions (n,3)`, `radii`, `landed`, etc. Droplets do
  not interact, so where they are emitted from is arbitrary per droplet.
- `Nozzle.emit(n, rng)` emits all droplets from a single `NozzleConfig.position`.
- `analysis.py` has `radial_distances` / `summarize` for a single spray spot.

---

## A — Deposition analysis (thickness + uniformity)

### A1 — Deposition / thickness map (the shared primitive)

**Goal.** Turn landing positions + droplet sizes into a 2-D **film thickness
field** on the target surface.

**Model.** Bin landed droplets onto a regular grid over the surface (`z = ground`).
For each cell, the deposited **wet volume** is the sum of droplet volumes
`(4/3)π r³` landing in it; the **wet thickness** is `volume / cell_area`. The
sprayed liquid is a **solution/dispersion**: only a fraction of its volume is the
solid material that stays after the solvent evaporates, so the **dry (cured)
thickness** = `wet · solids_fraction`. (No spreading/coalescence/run-off or
evaporation kinetics modelled — a droplet deposits its volume where it lands; call
this out.)

**Material change (prerequisite).** Add `solids_fraction: float = 1.0` to
`MaterialConfig` — the fraction of the sprayed solution that is solid material used
to prepare it (the rest is solvent that evaporates). Interpreted as a **volume
fraction** for the film-thickness model (document this; mass-fraction inputs can be
converted via densities later). `1.0` = pure liquid, wet == dry. Plumb it through
CLI (`--solids-fraction`), config (`SOLIDS_FRACTION`) and `.npz` storage
(back-compat default `1.0`), mirroring how `density`/`viscosity` are handled.

**API (new in `spraysim/analysis.py`).**
- `DepositionField` dataclass: `x_edges`, `y_edges`, `thickness` (2-D, metres, the
  **dry** thickness), `cell_size`, `cell_area`, plus helpers (`extent`, `mean`,
  `nonzero_mask`).
- `deposition_map(result, config, *, cell_size=None, extent=None) -> DepositionField`.
  Uses `config.material.solids_fraction` for the dry thickness. Defaults:
  `extent` = landing-point bounding box (padded); `cell_size` = a few × mean
  droplet spacing (or an explicit value). Uses
  `np.histogram2d(weights=volume)` for speed.

**Files.** `spraysim/config.py` (+ `MaterialConfig.solids_fraction`),
`spraysim/analysis.py` (+ exports), `run.py`, `main.sh`, `config/*.conf`,
`spraysim/storage.py`; `docs/material_properties.md`; `tests/`.

**Acceptance.**
- Total gridded volume equals `Σ (4/3)π r³` over landed droplets (to rounding).
- A known count of equal droplets over a known area gives the expected mean
  thickness analytically.

### A2 — Uniformity metrics

**Goal.** Quantify how even the coating is, computed on the A1 field (mass-based)
and optionally on per-cell **counts** (particle-based).

**Metrics (`UniformityStats` dataclass + `uniformity(field, *, roi=None,
coverage_threshold=...) -> UniformityStats`).**
- **CV** = std/mean of cell thickness over the ROI (lower = more uniform).
- **Christiansen Uniformity Coefficient** `CU = 1 − Σ|xᵢ−x̄| / (n·x̄)` (the standard
  coating/irrigation uniformity index; 1 = perfect).
- **Coverage fraction** = fraction of ROI cells with thickness ≥ threshold.
- **min/max, p10/p90, mean thickness.**
- ROI defaults to the field's non-zero region, or a user rectangle (needed for B,
  where edges of a raster are legitimately thinner).

**Files.** `spraysim/analysis.py` (+ export). Tests.

**Acceptance.**
- A synthetic *uniform* field → CV ≈ 0, CU ≈ 1, coverage ≈ 1.
- A single-spot spray (existing default run) → CV and CU reflect the radial
  fall-off (CU well below 1); metrics are stable across seeds within tolerance.

### A3 — Visualisation (optional but recommended)

- Add a **thickness heatmap** panel to `spraysim/plots.py` (or a dedicated
  `plot_deposition`), and a thickness/uniformity page to `analysis/report.py`.
- Print CV / CU / coverage in the `run.py` stats block and store them in the
  `.npz` (so saved runs carry their uniformity).

**Acceptance.** Report renders a heatmap; stats block shows CU/CV/coverage.

---

## B — G-code path spraying

**Goal.** Spray while the nozzle travels a path from a `.gcode` file, depositing a
2-D pattern (line, raster, arbitrary toolpath) rather than a single spot.

**Key enabling fact.** Droplets are independent, so a moving source is just a
batch of droplets emitted from many positions/times. We can generate the whole
path's droplets and feed them to the **existing** integrator in one batch — no
change to the integrator core.

### B1 — G-code parser (`spraysim/gcode.py`)

Parse a pragmatic subset into motion segments:
- `Move(start, end, feed, spray_on)` with `start/end` in metres (surface frame),
  `feed` in m/s.
- Supported: `G0` (rapid → **spray off / travel**), `G1` (linear → **spray on**),
  `X Y Z` coordinates, `F` feed rate, `G90/G91` (abs/rel), `G20/G21`
  (inch/mm), comments (`;`, `( )`). Arcs `G2/G3` deferred (approximate as line, or
  raise).
- **Spray on/off:** `G1` = spray on, `G0` = travel (spray off). Fixed convention.
- **Units:** coordinates are **mm** by default (`G21`); converted to metres
  internally. `G20` (inch) supported for completeness.
- Helpers: `total_spray_time(moves)`, `bounds(moves)`.

**Files.** new `spraysim/gcode.py`; unit tests with sample programs.

**Acceptance.** Parses a hand-written raster; abs/rel and mm/inch handled; spray
on/off segmentation correct; total spray length/time match by hand.

### B2 — Path emitter (generalise `Nozzle`)

- Refactor `Nozzle.emit` to `emit_from(position, n, rng, *,
  carriage_velocity=None)` — emit `n` droplets from an arbitrary position, adding
  the nozzle's **carriage velocity** to each droplet's launch velocity (a moving
  nozzle throws droplets along the travel direction; important at high feed).
- `emit_path(moves, rng) -> (positions, velocities, radii)`: for each spray-on
  segment, droplets = `flow_rate · segment_time / V_droplet`; place them at
  positions interpolated along the segment (uniform in path length), each carrying
  that segment's carriage velocity and the current `Z` as standoff height.
- Total count = `flow_rate · total_spray_time / V_droplet`, capped by
  `max_droplets` (same balance as today, `spray_duration → path spray time`).

**Files.** `spraysim/nozzle.py`. Tests.

**Acceptance.** A 2-segment line deposits droplets along that line; disabling
carriage velocity vs enabling shifts landings downstream by ~`v_carriage · t_fall`.

### B3 — Config & CLI

- New `PathConfig` dataclass: `gcode` (file path or inline text), `default_feed`,
  `include_carriage_velocity=True`, spray-signal convention. `SimConfig.path:
  PathConfig | None = None` — when `None`, current single-spot behaviour; when set,
  path spraying (and `spray_duration` is derived from the path).
- CLI: `--gcode FILE` (+ `--feed`, `--no-carriage-velocity`); `main.sh` key
  `GCODE=` / `FEED=`; a `raster.gcode` example under `config/` or `examples/`.

**Files.** `spraysim/config.py`, `run.py`, `main.sh`, an example `.gcode`.

**Acceptance.** `./main.sh` with a `GCODE=` set runs a path; without it, unchanged.

### B4 — Simulator wiring & persistence

- In `Simulator.run()`: if `cfg.path` is set, build the droplet batch via
  `emit_path` (deriving the count from path spray time); else use `emit` as now.
- `SimResult`/`storage`: record the path (file name + parsed segments or bounds)
  and per-droplet source position (optional) so reports can draw the toolpath.
- `analysis/report.py`: overlay the toolpath on the landing/thickness plots.

**Files.** `spraysim/simulator.py`, `spraysim/storage.py`, `analysis/report.py`.

**Acceptance.** A raster `.gcode` produces a filled rectangle; `.npz` reloads with
the path; report shows toolpath + thickness map.

### B5 — Validation (ties A and B together)

- `analysis/validate.py` + tests: a **dense raster** with adequate overlap yields
  a **near-uniform** interior film (high CU, low CV in the ROI); a **sparse**
  raster shows periodic banding (lower CU) — the metrics detect the difference.

---

## Risks & open questions

- **Droplet-count explosion.** Long paths × derived flow can blow past
  `max_droplets`. Mitigate with the existing cap plus, if needed, a per-segment
  density / statistical down-sampling with volume re-weighting (keep total mass).
- **G-code dialects.** Many exist. Ship a **documented supported subset** and
  raise clearly on unsupported codes (esp. `G2/G3` arcs — decide: linearise or
  reject initially).
- **Moving-nozzle aerodynamics.** We add the *kinematic* carriage velocity only;
  air entrainment / wake from motion is not modelled. State this.
- **Deposition physics.** Thickness assumes deposit-where-it-lands (no spreading,
  coalescence, run-off, or evaporation). Fine for coating-uniformity comparisons;
  document the limitation.

**Resolved decisions.**
1. Spray on/off — **`G1` = on, `G0` = travel** (no `M3/M5` needed).
2. G-code coordinate units — **mm** (`G21`) by default.
3. `solids_fraction` — **a `MaterialConfig` property** (volume fraction of solids
   in the prepared solution; drives dry film thickness).

---

## Suggested sequencing

1. **A1** deposition/thickness map → **A2** uniformity metrics (independent,
   immediately useful on current single-spot runs, and the acceptance tool for B).
2. **A3** visualisation/reporting hooks.
3. **B1** G-code parser → **B2** path emitter → **B3** config/CLI → **B4** wiring
   & persistence.
4. **B5** end-to-end validation (raster uniformity), closing the loop with A.

Each step: implement → `pytest` + `analysis/validate.py` → commit file-by-file
with descriptive messages (per the project convention).
