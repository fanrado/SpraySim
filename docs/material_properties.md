# Material Properties

This is the second of the **two inputs that dictate a spray result** (the other
is the [sprayer](sprayer_parameters.md)). It describes the *sprayed liquid* — the
substance leaving the nozzle. The simulation is **not limited to water**: any
liquid can be sprayed, and its properties feed directly into the physics.

The single material property the model currently uses is **density**
(`ρ_liquid`, kg/m³). It is defined in `spraysim/config.py` as `MaterialConfig`
and backed by a small registry in `spraysim/materials.py`.

---

## Configuring the material

Pick a liquid **by name** (density comes from the registry) or give **any name
plus an explicit density** for something not listed.

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Material name | `MATERIAL` | `--material` | `MaterialConfig.name` | name | `water` | Selects a registry density and labels the run/output. |
| Density override | `DENSITY` | `--density` | `MaterialConfig.density` | kg/m³ | *(empty → registry default)* | Overrides the density for a custom or off-registry liquid. |

```bash
# By name (registry density)
python run.py --material diesel

# Custom liquid: any name + explicit density
python run.py --material glycol --density 1113
```

```python
from spraysim import MaterialConfig, SimConfig, Simulator
material = MaterialConfig(name="ethanol", density=789.0)
result = Simulator(SimConfig(material=material)).run()
```

An unknown material name **without** a `--density` override fails with a clear
error listing the known materials.

---

## Built-in material registry

Nominal densities near 15–20 °C, from `spraysim/materials.py`.

| Material | Density (kg/m³) | Notes |
|----------|-----------------|-------|
| `water` | 1000 | Reference liquid *(default)* |
| `seawater` | 1025 | Slightly denser than fresh water |
| `ethanol` | 789 | Light alcohol |
| `methanol` | 792 | Light alcohol |
| `acetone` | 784 | Light solvent |
| `gasoline` | 745 | Lightest listed fuel |
| `kerosene` | 810 | Fuel |
| `diesel` | 832 | Fuel |
| `olive_oil` | 915 | Viscous oil (only density is modelled) |
| `glycerin` | 1260 | Densest listed liquid |

For anything else, pass `--density` (or set `MaterialConfig.density`) with the
value for your liquid and temperature.

---

## How density shapes the spray

Density enters the physics in **two independent places**, so changing the liquid
changes the result even with an identical sprayer.

### 1. Exit hydraulics (Torricelli)

```
v_ideal = √(2 · ΔP / ρ_liquid)
```

A denser liquid is harder to accelerate through the orifice, so for the **same
pressure**:

| Quantity | Scaling | Effect of denser liquid |
|----------|---------|--------------------------|
| Exit speed `v` | `∝ 1/√ρ_liquid` | Slower |
| Flow rate `Q` | `∝ 1/√ρ_liquid` | Lower |
| Droplet count `N` | `∝ 1/√ρ_liquid` | Fewer (less volume flows in the same time) |

Example at 3 bar through a `full_cone` tip: water (1000 kg/m³) leaves at
~23.3 m/s, ethanol (789 kg/m³) at ~26.2 m/s (`v ∝ 1/√ρ`).

### 2. Droplet mass / aerodynamic drag

Each droplet is a sphere of the liquid, so mass `m = ρ_liquid · (4/3)π r³`. The
per-droplet drag factor is:

```
a_drag = -k · |v| · v,     k = 3 · ρ_air · C_d / (8 · ρ_liquid · r)
```

A denser liquid gives a **smaller `k`** (heavier droplet, more inertia per unit
frontal area), so it decelerates *less* and, once airborne, tends to travel
*further* before landing.

### Net effect

The two effects pull in opposite directions on range — a denser liquid leaves
**slower** (shorter throw) but resists drag **better** (longer throw). The
simulation resolves the balance for you; there is no simple closed form for the
combined landing distance, which is exactly why running the sim is useful.

---

## Scope and limitations

- **Only density is modelled.** Viscosity and surface tension are *not* part of
  the current physics. They would matter for real atomisation — how the sheet
  breaks up and what droplet sizes actually form — but here the droplet size
  distribution is an **input** (see the sprayer parameters), not derived from the
  liquid. Listing a viscous liquid like `olive_oil` only changes its density.
- Density is treated as constant (incompressible, isothermal). No temperature or
  evaporation effects.
- The material name and density are saved into the `.npz` output, so every run is
  fully reproducible and self-describing.

See [sprayer_parameters.md](sprayer_parameters.md) for the other half of the
inputs.
