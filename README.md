# SpraySim

A physics-based **spray / droplet particle simulator**. A nozzle emits droplets
into a velocity cone; each droplet is integrated under **gravity** and
**quadratic aerodynamic drag** until it lands on the ground plane. The output is
numerical statistics plus static summary plots (no animation).

## Physics model

Each droplet is a sphere of water. Two forces act on it:

- **Gravity:** `a = -g ẑ`
- **Aerodynamic drag:** `F = -½ ρ_air C_d A |v| v`, with cross-section
  `A = π r²` and mass `m = ρ_water · (4/3) π r³`.

Combining drag into a per-droplet factor gives the acceleration used by the
integrator:

```
a_drag = -k |v| v,   k = 3 ρ_air C_d / (8 ρ_water r)
```

Smaller droplets have a larger `k`, so they decelerate faster and travel less —
the physical reason a fine mist stays close while big drops fly further.

Integration uses **semi-implicit (symplectic) Euler** with a fixed timestep, and
the ground impact is found by linearly interpolating the crossing point within
the final step.

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

| Key         | Meaning                          | `run.py` flag  |
|-------------|----------------------------------|----------------|
| `DROPLETS`  | number of droplets               | `--droplets`   |
| `SPEED`     | nozzle exit speed (m/s)          | `--speed`      |
| `CONE`      | cone half-angle (degrees)        | `--cone`       |
| `HEIGHT`    | nozzle height (m)                | `--height`     |
| `RADIUS_MM` | mean droplet radius (mm)         | `--radius-mm`  |
| `DT`        | integration timestep (s)         | `--dt`         |
| `SEED`      | RNG seed                         | `--seed`       |
| `OUT`       | output figure path               | `--out`        |
| `NO_PLOT`   | `true` = stats only, no figure   | `--no-plot`    |

Shipped presets: `default`, `fine_mist` (small slow droplets, wide cone) and
`big_drops` (large fast droplets, narrow cone).

### Calling `run.py` directly

`main.sh` just translates a config into a `run.py` invocation, so you can also
run the Python CLI directly for ad-hoc experiments:

```bash
python run.py                              # default spray -> output/spray_summary.png
python run.py --droplets 8000 --speed 12 --cone 35
python run.py --height 2.0 --radius-mm 0.6 --out output/big_drops.png
python run.py --no-plot                    # print statistics only
```

Either way, the run prints a JSON block of statistics (coverage radius, mean
flight time, impact speed, …) and writes a 2×2 figure:

1. **Trajectories (side view)** — sampled droplet paths, coloured by radius.
2. **Landing pattern (top view)** — where droplets land, coloured by impact speed.
3. **Radial deposition profile** — histogram of distance from the spray axis.
4. **Droplet size distribution** — the sampled radii.

## Use as a library

```python
from spraysim import SimConfig, NozzleConfig, Simulator, analysis, plots

config = SimConfig(n_droplets=5000, nozzle=NozzleConfig(exit_speed=10.0))
result = Simulator(config).run()
stats = analysis.summarize(result, config)
print(stats.as_dict())
plots.save_figure(result, config, "output/custom.png")
```

## Project layout

```
config/          # *.conf presets (KEY=value) — the inputs you edit
main.sh          # launcher: loads a config and runs the simulation
spraysim/
  config.py      # dataclasses: PhysicsConfig, NozzleConfig, SimConfig
  nozzle.py      # samples initial positions, cone velocities, droplet radii
  simulator.py   # vectorised integrator -> SimResult
  analysis.py    # derived statistics (coverage, flight time, ...)
  plots.py       # static matplotlib summary figure
run.py           # Python CLI entry point (called by main.sh)
tests/           # pytest sanity checks (incl. vacuum free-fall vs analytic)
```

## Test

```bash
pytest
```
