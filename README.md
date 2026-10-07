# SpraySim

A physics-based **spray / droplet particle simulator**. A nozzle emits droplets
into a velocity cone; each droplet is integrated under **gravity** and
**Reynolds-dependent aerodynamic drag** until it lands on the ground plane. The
output is numerical statistics, static summary plots, and a self-describing
`.npz` data archive per run.

Droplet count and exit speed are **derived** from the nozzle pressure, orifice
and shape (not typed in), the sprayed liquid is a configurable **material**, and
the droplet-size distribution is selectable.

## Documentation

- [**Physics model**](docs/physics.md) — forces, drag, hydraulics, size
  distribution, and how the material's density enters.
- [**Sprayer parameters**](docs/sprayer_parameters.md) — every nozzle / run input
  (geometry, pressure, orifice, shape, distribution, drag model, controls).
- [**Material properties**](docs/material_properties.md) — the sprayed liquid, its
  density and viscosity, and the built-in registry.
- [**Analysis**](analysis/README.md) — turn saved `.npz` runs into a PDF report.

## Environment setup

Requires **Python 3.9+** with NumPy and Matplotlib (only third-party deps).

```bash
python -m venv .venv && source .venv/bin/activate   # optional
pip install -r requirements.txt
```

## Run

```bash
./main.sh
```

That runs the simulation with `config/default.conf` and writes the results to
`output/tmp/`. To use another config or keep the results, name them:

```bash
./main.sh fine_mist --run-name my_run     # config/fine_mist.conf -> output/my_run/
./main.sh --list                          # list the available configs
```

A config is a plain `KEY=value` file in `config/`; copy one to make a new
preset. Any `run.py` flag after the config name is passed straight through and
overrides the config (`./main.sh default --cone 15 --no-plot`).

### Parameters

Every parameter `main.sh` reads from the config, with its `run.py` flag. The
physics behind each one is in [sprayer_parameters.md](docs/sprayer_parameters.md)
and [material_properties.md](docs/material_properties.md).

| Config key | `run.py` flag | Default | Description |
|------------|---------------|---------|-------------|
| **Material (sprayed liquid)** | | | |
| `MATERIAL` | `--material` | `water` | Liquid, by name: `water`, `seawater`, `ethanol`, `methanol`, `acetone`, `toluene`, `gasoline`, `kerosene`, `diesel`, `olive_oil`, `glycerin`, or any name with `DENSITY`. |
| `DENSITY` | `--density` | *(registry)* | Liquid density override, kg/m³. |
| `VISCOSITY` | `--viscosity` | *(registry)* | Liquid dynamic viscosity override, Pa·s (reported, not yet used in flight). |
| `SOLIDS_FRACTION` | `--solids-fraction` | `1.0` | Volume fraction of solids in the solution; dry film thickness = wet × this. |
| **Hydraulics (set exit speed, flow and droplet count)** | | | |
| `PRESSURE_BAR` | `--pressure-bar` | `3.0` | Nozzle pressure, bar. |
| `ORIFICE_MM` | `--orifice-mm` | `0.8` | Orifice diameter, mm. |
| `NOZZLE_SHAPE` | `--shape` | `full_cone` | `sharp_orifice`, `rounded_orifice`, `full_cone`, `hollow_cone` or `flat_fan`. |
| `SPRAY_DURATION` | `--spray-duration` | `0.15` | Seconds the nozzle is open (fixed-spot runs); scales the droplet count. |
| `DROPLETS` | `--droplets` | *(empty → derived)* | Pin the droplet count instead of deriving it. |
| **Droplet size distribution** | | | |
| `DISTRIBUTION` | `--distribution` | `lognormal` | `normal` or `lognormal`. |
| `MEAN_RADIUS_MM` | `--mean-radius-mm` | `0.4` | Mean droplet radius, mm. |
| `RADIUS_STD_MM` | `--radius-std-mm` | `0.12` | Standard deviation of the radius, mm. |
| **Toolpath (optional: spray while moving)** | | | |
| `GCODE` | `--gcode` | *(empty → fixed spot)* | G-code file; `G1` = spray, `G0` = travel. See [Toolpaths](#toolpaths-g-code). |
| `FEED` | `--feed` | *(empty → program `F`)* | Feed-rate override, mm/min. |
| `STANDOFF_MM` | `--standoff-mm` | `150` | Nozzle height above the surface, mm (path runs). |
| — | `--no-carriage-velocity` | off | Do not add the nozzle travel velocity to the droplets. |
| **Physics** | | | |
| `DRAG_MODEL` | `--drag-model` | `clift_gauvin` | `clift_gauvin` (Reynolds-dependent drag) or `constant` (fixed C_d). |
| **Geometry / integration** | | | |
| `CONE` | `--cone` | `25.0` | Spray-cone half-angle, degrees. |
| `HEIGHT` | `--height` | `1.5` | Nozzle height, m (fixed-spot runs). |
| `SPEED_SPREAD` | `--speed-spread` | `0.15` | Relative spread of droplet speed about the exit speed. |
| `DT` | `--dt` | `0.001` | Integration timestep, s. |
| `SEED` | `--seed` | `42` | RNG seed. |
| **Output** | | | |
| `RUN_NAME` | `--run-name` | *(empty → `tmp`)* | Folder for this run's files: `output/<RUN_NAME>/`. |
| `OUT` | `--out` | `spray_summary.png` | Summary figure, file name inside the run folder. |
| `NO_PLOT` | `--no-plot` | `false` | Skip the figure. |
| `DATA` | `--data` | `spray_data.npz` | Result archive, file name inside the run folder. |
| `NO_DATA` | `--no-data` | `false` | Skip the archive. |

`main.sh` only turns the config into a `run.py` call, so `python run.py --help`
lists the same flags and `python run.py --material diesel --pressure-bar 5` runs
without a config.

### Outputs

Each run writes into its own folder `output/<RUN_NAME>/` (`output/tmp/` when
unnamed, overwritten by the next unnamed run), prints a JSON block of statistics
— including **deposition & uniformity** (dry film thickness, CV, Christiansen
CU, coverage) — and writes:

- a **2×2 summary figure** (`--out`, skip with `--no-plot`);
- a compressed **`.npz` archive** (`--data`, skip with `--no-data`) holding the
  per-droplet arrays and the full config, so a run reloads exactly:

  ```python
  from spraysim import storage
  result, config = storage.load_result("output/tmp/spray_data.npz")
  ```

- optionally a multi-page **PDF report** across saved runs (all run folders):
  `python analysis/report.py`;
- a fit of a measured wavelength-shifter efficiency vs. deposited density with
  the simulated coating as the model (`python analysis/fit_ptp_efficiency.py`,
  see [analysis/README.md](analysis/README.md)).

## Toolpaths (G-code)

Set `GCODE=` to spray along a path instead of from a fixed spot; the deposited
film and its uniformity (CV, Christiansen CU, coverage) are reported per run.
Two helpers make and check toolpaths:

```bash
python svg_to_gcode.py drawing.svg --fit-box-mm 0 0 120 120 --closed-loop   # SVG -> G-code
./run_svg_to_gcode.sh drawing.svg                                           # same, preset box
python gcode_to_svg.py path.gcode --show-travel                             # G-code -> SVG (check)
```

`svg_to_gcode.py` linearises the SVG `<path>` elements (lines, Beziers, arcs).
`--fit-box-mm` fits and centres the artwork in a mm box; `--closed-loop` appends
a spray-off return pass so repeated passes start where they ended
(`--return-feed` sets its speed); `--home` and `--z-offset-mm` optionally add a
homing move and a Z height (otherwise `STANDOFF_MM` sets the height). See
`--help` for scale, units, tolerance and Y-flip options. `gcode_to_svg.py` is the
exact inverse and round-trips the artwork, so it is the validator.

## Use as a library

```python
from spraysim import SimConfig, NozzleConfig, MaterialConfig, Simulator, analysis

nozzle = NozzleConfig(pressure=5.0e5, orifice_diameter=1.0e-3, shape="flat_fan",
                      distribution="normal", mean_radius=3.0e-4, radius_std=8.0e-5)
material = MaterialConfig(name="diesel", density=832.0)
config = SimConfig(nozzle=nozzle, material=material, spray_duration=0.2)
result = Simulator(config).run()
print(analysis.summarize(result, config).as_dict())
```

## Project layout

```
docs/            # physics.md + input reference docs
analysis/        # offline analysis of saved runs -> PDF report + validation
config/          # *.conf presets (KEY=value) — the inputs you edit
examples/        # example G-code toolpaths (raster.gcode, raster_120mm_pitch12.gcode)
output/          # one folder per run: output/<RUN_NAME>/ (git-ignored)
main.sh          # launcher: loads a config and runs the simulation
run.py           # Python CLI entry point (called by main.sh)
svg_to_gcode.py  # SVG artwork -> G-code toolpath; gcode_to_svg.py is its inverse
spraysim/        # the package (config, hydraulics, drag, gcode, nozzle, simulator, ...)
tests/           # pytest sanity + physics-validation checks
```

## Platform

Developed and tested on **macOS 15 (Darwin, Apple Silicon)** with **CPython
3.13**, NumPy 2.3 and Matplotlib 3.10. It is pure Python plus NumPy/Matplotlib
with no OS-specific code, so it is expected to run on Linux and Windows as well
(the `main.sh` launcher needs a POSIX shell; on Windows use `run.py` directly).

## Test

```bash
pytest
```
