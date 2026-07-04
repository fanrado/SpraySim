# SpraySim

A physics-based **spray / droplet particle simulator**. A nozzle emits droplets
into a velocity cone; each droplet is integrated under **gravity** and
**quadratic aerodynamic drag** until it lands on the ground plane. The output is
numerical statistics plus static summary plots (no animation).

A spray result is dictated by **two inputs**, each documented in detail:

- [**Sprayer parameters**](docs/sprayer_parameters.md) — geometry, pressure,
  orifice, nozzle shape, droplet-size distribution and run controls.
- [**Material properties**](docs/material_properties.md) — the sprayed liquid and
  its density.

## Physics model

Each droplet is a sphere of the sprayed liquid. Two forces act on it:

- **Gravity:** `a = -g ẑ`
- **Aerodynamic drag:** `F = -½ ρ_air C_d A |v| v`, with cross-section
  `A = π r²` and mass `m = ρ_liquid · (4/3) π r³`.

Combining drag into a per-droplet factor gives the acceleration used by the
integrator:

```
a_drag = -k |v| v,   k = 3 ρ_air C_d / (8 ρ_liquid r)
```

Smaller droplets have a larger `k`, so they decelerate faster and travel less —
the physical reason a fine mist stays close while big drops fly further.

Integration uses **semi-implicit (symplectic) Euler** with a fixed timestep, and
the ground impact is found by linearly interpolating the crossing point within
the final step.

### Material (sprayed liquid)

The simulation is not limited to water: the sprayed liquid is a configurable
**material** whose density `ρ_liquid` feeds both the droplet mass in the drag
term above and the Torricelli exit speed / flow rate below. Pick a liquid by
name from a small registry (`water`, `seawater`, `ethanol`, `methanol`,
`acetone`, `gasoline`, `kerosene`, `diesel`, `olive_oil`, `glycerin`) or give
any name together with an explicit density. A denser liquid exits more slowly
for the same pressure (`v ∝ 1/√ρ_liquid`) but carries more momentum per drop.

Each material also carries a **dynamic viscosity** (`μ_liquid`, Pa·s), defaulting
to water's value. It is stored and reported but not yet coupled to the flight
physics — it is the hook for future viscous/atomisation effects. See
[docs/material_properties.md](docs/material_properties.md).

### Nozzle hydraulics (droplet count & exit speed)

Rather than typing a droplet count and speed directly, both are **derived** from
the nozzle pressure, orifice size and shape via incompressible orifice flow
(Bernoulli / Torricelli):

```
ideal velocity   v_ideal = √(2 ΔP / ρ)
exit speed       v       = C_v · v_ideal
volumetric flow  Q       = C_d · A · v_ideal          (A = orifice area)
droplet count    N       = Q · spray_duration / V_droplet
```

`C_d` (discharge) and `C_v` (velocity) coefficients come from the nozzle
**shape** (`sharp_orifice`, `rounded_orifice`, `full_cone`, `hollow_cone`,
`flat_fan`). `V_droplet = (4/3)π·E[r³]` is the mean droplet volume, so for a
fixed flow, finer droplets mean *many* more of them (`N ∝ 1/r³`).

### Droplet size distribution

Droplet radii are drawn from a configurable **`normal`** (Gaussian) or
**`lognormal`** distribution, both parameterised by a linear-space mean radius
and standard deviation. `E[r³]` (which sets the droplet count) is computed
analytically from these for each distribution.

## Install

```bash
pip install -r requirements.txt
```

## Run

The recommended way to run is through **`main.sh`**, which loads parameters from
a config file so you never have to type flags:

```bash
./main.sh                        # uses config/default.conf
./main.sh fine_mist              # uses config/fine_mist.conf (by name)
./main.sh config/big_drops.conf  # explicit path also works
./main.sh --list                 # list available configs
./main.sh default --no-plot      # extra flags pass straight through to run.py
```

`main.sh` resolves the project root itself, so it works from any directory.

### Configs

Inputs live in `config/*.conf` — plain shell `KEY=value` files. Edit these
instead of passing command-line flags. To make a new preset, copy one:

```bash
cp config/default.conf config/my_run.conf
# edit my_run.conf, then:
./main.sh my_run
```

| Key              | Meaning                                        | `run.py` flag       |
|------------------|------------------------------------------------|---------------------|
| `MATERIAL`       | sprayed liquid (name from the registry)        | `--material`        |
| `DENSITY`        | *optional* liquid density override (kg/m³)     | `--density`         |
| `VISCOSITY`      | *optional* dynamic viscosity override (Pa·s)   | `--viscosity`       |
| `PRESSURE_BAR`   | nozzle pressure (bar)                          | `--pressure-bar`    |
| `ORIFICE_MM`     | orifice diameter (mm)                          | `--orifice-mm`      |
| `NOZZLE_SHAPE`   | shape → discharge/velocity coefficients        | `--shape`           |
| `SPRAY_DURATION` | seconds the nozzle is open (scales count)      | `--spray-duration`  |
| `DISTRIBUTION`   | `normal` or `lognormal`                        | `--distribution`    |
| `MEAN_RADIUS_MM` | mean droplet radius (mm)                       | `--mean-radius-mm`  |
| `RADIUS_STD_MM`  | std of droplet radius (mm)                     | `--radius-std-mm`   |
| `CONE`           | cone half-angle (degrees)                      | `--cone`            |
| `HEIGHT`         | nozzle height (m)                              | `--height`          |
| `SPEED_SPREAD`   | relative std of speed about the exit speed     | `--speed-spread`    |
| `DT`             | integration timestep (s)                       | `--dt`              |
| `SEED`           | RNG seed                                        | `--seed`            |
| `OUT`            | output figure path                             | `--out`             |
| `NO_PLOT`        | `true` = stats only, no figure                 | `--no-plot`         |
| `DATA`           | output `.npz` data path (arrays + config)      | `--data`            |
| `NO_DATA`        | `true` = skip saving the `.npz` data           | `--no-data`         |
| `DROPLETS`       | *optional* explicit count (empty = derive)     | `--droplets`        |

Droplet **count** and **exit speed** are derived from `PRESSURE_BAR`,
`ORIFICE_MM` and `NOZZLE_SHAPE` (see *Nozzle hydraulics* above); leave
`DROPLETS` empty to use that derivation, or set it to pin an explicit count.

Shipped presets: `default` (water), `fine_mist` (low-pressure hollow-cone
atomiser spraying ethanol, small Gaussian droplets) and `big_drops`
(high-pressure flat-fan water, large lognormal droplets).

### Calling `run.py` directly

`main.sh` just translates a config into a `run.py` invocation, so you can also
run the Python CLI directly for ad-hoc experiments:

```bash
python run.py                                              # default nozzle (water)
python run.py --material diesel --pressure-bar 5 --shape flat_fan
python run.py --material glycol --density 1113            # custom liquid by density
python run.py --distribution normal --mean-radius-mm 0.3 --radius-std-mm 0.08
python run.py --droplets 5000 --no-plot                   # pin count, stats only
```

Either way, the run prints a JSON block of statistics (coverage radius, mean
flight time, impact speed, …), writes an `.npz` data file (see *Output data*
below) and writes a 2×2 figure:

1. **Trajectories (side view)** — sampled droplet paths, coloured by radius.
2. **Landing pattern (top view)** — where droplets land, coloured by impact speed.
3. **Radial deposition profile** — histogram of distance from the spray axis.
4. **Droplet size distribution** — the sampled radii.

## Output data (`.npz`)

Every run also saves the full result to a compressed NumPy archive (`.npz`, a
dictionary of named arrays) so you can re-analyse a run without recomputing it.
The archive holds the per-droplet arrays (`landing_positions`, `flight_times`,
`impact_speeds`, `radii`, `launch_speeds`, `landed`), the sampled trajectories
(flattened into `traj_points` + `traj_lengths`), the derived hydraulics, and the
complete run config (`cfg_*` keys) so it is fully self-describing.

Load it back with the raw NumPy API for ad-hoc analysis, or with the helper that
reconstructs the `SimResult` and `SimConfig`:

```python
import numpy as np
z = np.load("output/spray_data.npz")
print(z["radii"].mean(), z["landing_positions"].shape)

from spraysim import storage
result, config = storage.load_result("output/spray_data.npz")
# `config` reproduces the run exactly: Simulator(config).run() gives the same result.
```

Pass `--no-data` (or `NO_DATA=true`) to skip writing it.

## Use as a library

```python
from spraysim import (SimConfig, NozzleConfig, MaterialConfig,
                      Simulator, analysis, plots, storage)

nozzle = NozzleConfig(pressure=5.0e5, orifice_diameter=1.0e-3, shape="flat_fan",
                      distribution="normal", mean_radius=3.0e-4, radius_std=8.0e-5)
material = MaterialConfig(name="diesel", density=832.0)  # spray something other than water
config = SimConfig(nozzle=nozzle, material=material, spray_duration=0.2)  # count derived
result = Simulator(config).run()

print("exit speed:", result.exit_speed, "m/s  droplets:", result.n)
print(analysis.summarize(result, config).as_dict())
plots.save_figure(result, config, "output/custom.png")
storage.save_result(result, config, "output/custom.npz")  # reload later, no recompute
```

## Project layout

```
docs/            # reference docs: sprayer_parameters.md, material_properties.md
config/          # *.conf presets (KEY=value) — the inputs you edit
main.sh          # launcher: loads a config and runs the simulation
spraysim/
  config.py      # dataclasses: PhysicsConfig, MaterialConfig, NozzleConfig, SimConfig
  materials.py   # sprayed-liquid registry (name -> density)
  hydraulics.py  # pressure/orifice/shape -> exit speed, flow rate, droplet count
  nozzle.py      # samples initial positions, cone velocities, droplet radii
  simulator.py   # vectorised integrator -> SimResult
  analysis.py    # derived statistics (coverage, flight time, ...)
  plots.py       # static matplotlib summary figure
  storage.py     # save/load a run to a .npz archive (arrays + config)
run.py           # Python CLI entry point (called by main.sh)
tests/           # pytest sanity checks (incl. vacuum free-fall vs analytic)
```

## Test

```bash
pytest
```
