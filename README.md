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

The recommended way is through **`main.sh`**, which loads parameters from a config
file so you never have to type flags:

```bash
./main.sh                        # uses config/default.conf
./main.sh fine_mist              # uses config/fine_mist.conf (by name)
./main.sh config/big_drops.conf  # explicit path also works
./main.sh --list                 # list available configs
./main.sh default --no-plot      # extra flags pass straight through to run.py
```

Inputs live in `config/*.conf` — plain shell `KEY=value` files; copy one to make a
new preset. Every key is documented in
[sprayer_parameters.md](docs/sprayer_parameters.md) and
[material_properties.md](docs/material_properties.md). Shipped presets: `default`
(water), `fine_mist` (ethanol, small droplets), `big_drops` (large droplets), and
`raster` (moves the nozzle along a **G-code toolpath** to build a uniform
coating).

Set `GCODE=` (or `--gcode file.gcode`) to spray while moving along a path
(`G1` = spray, `G0` = travel) instead of a fixed spot; the deposited film and its
uniformity (CV / Christiansen CU / coverage) are reported per run.

`main.sh` just translates a config into a `run.py` invocation, so you can also run
the CLI directly:

```bash
python run.py                                    # default nozzle (water)
python run.py --material diesel --pressure-bar 5 --shape flat_fan
python run.py --distribution normal --mean-radius-mm 0.3 --radius-std-mm 0.08
python run.py --help                             # all flags
```

### Outputs

Each run prints a JSON block of statistics — including **deposition & uniformity**
(dry film thickness, CV, Christiansen CU, coverage) — and writes:

- a **2×2 summary figure** (`--out`, skip with `--no-plot`);
- a compressed **`.npz` archive** (`--data`, skip with `--no-data`) holding the
  per-droplet arrays and the full config, so a run reloads exactly:

  ```python
  from spraysim import storage
  result, config = storage.load_result("output/spray_data.npz")
  ```

- optionally a multi-page **PDF report** across saved runs:
  `python analysis/report.py`.

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
examples/        # example G-code toolpaths (e.g. raster.gcode)
main.sh          # launcher: loads a config and runs the simulation
spraysim/        # the package (config, hydraulics, drag, gcode, nozzle, simulator, ...)
run.py           # Python CLI entry point (called by main.sh)
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
