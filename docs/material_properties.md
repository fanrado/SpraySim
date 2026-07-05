# Material Properties

This is the second of the **two inputs that dictate a spray result** (the other
is the [sprayer](sprayer_parameters.md)). It describes the *sprayed liquid* — the
substance leaving the nozzle. The simulation is **not limited to water**: any
liquid can be sprayed, and its properties feed directly into the physics.

A material carries these properties, defined in `spraysim/config.py` as
`MaterialConfig` and backed by a registry in `spraysim/materials.py`:

- **density** (`ρ_liquid`, kg/m³) — the property that currently drives the
  physics (droplet mass and exit hydraulics).
- **dynamic viscosity** (`μ_liquid`, Pa·s) — carried for custom-liquid
  definitions and reporting; defaults to water's value and is **not yet coupled**
  to the flight physics (see *Scope and limitations*).
- **solids fraction** — the volume fraction of solid material in the prepared
  **solution** (the rest is solvent that evaporates). It does not affect flight;
  it sets the **dry deposited film thickness** = wet thickness × solids_fraction
  (see the deposition analysis). `1.0` = pure liquid (wet == dry).

---

## Configuring the material

Pick a liquid **by name** (density comes from the registry) or give **any name
plus an explicit density** for something not listed.

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Material name | `MATERIAL` | `--material` | `MaterialConfig.name` | name | `water` | Selects registry properties and labels the run/output. |
| Density override | `DENSITY` | `--density` | `MaterialConfig.density` | kg/m³ | *(empty → registry default)* | Overrides the density for a custom or off-registry liquid. |
| Viscosity override | `VISCOSITY` | `--viscosity` | `MaterialConfig.viscosity` | Pa·s | *(empty → registry value; custom → water)* | Overrides the dynamic viscosity. |
| Solids fraction | `SOLIDS_FRACTION` | `--solids-fraction` | `MaterialConfig.solids_fraction` | volume fraction | `1.0` | Fraction of the solution that is solids; sets dry film thickness. |

```bash
# By name (registry density + viscosity)
python run.py --material diesel

# Custom liquid: any name + explicit density (viscosity defaults to water's)
python run.py --material glycol --density 1113

# Custom liquid with its own viscosity too
python run.py --material glycol --density 1113 --viscosity 0.0161
```

```python
from spraysim import MaterialConfig, SimConfig, Simulator
material = MaterialConfig(name="ethanol", density=789.0, viscosity=1.2e-3)
result = Simulator(SimConfig(material=material)).run()
```

An unknown material name **without** a `--density` override fails with a clear
error listing the known materials. Viscosity, by contrast, is optional for any
liquid: if you do not give one, a registered material uses its registry value and
a custom liquid falls back to **water's viscosity**.

---

## Built-in material registry

Nominal values near ~20 °C, from `spraysim/materials.py`. Viscosity is dynamic
(Pa·s); `1 mPa·s = 1e-3 Pa·s`.

| Material | Density (kg/m³) | Viscosity (Pa·s) | Notes |
|----------|-----------------|------------------|-------|
| `water` | 1000 | 1.00e-3 | Reference liquid *(default)* |
| `seawater` | 1025 | 1.08e-3 | Slightly denser than fresh water |
| `ethanol` | 789 | 1.20e-3 | Light alcohol |
| `methanol` | 792 | 0.59e-3 | Light alcohol |
| `acetone` | 784 | 0.32e-3 | Light, low-viscosity solvent |
| `gasoline` | 745 | 0.60e-3 | Lightest listed fuel |
| `kerosene` | 810 | 1.50e-3 | Fuel |
| `diesel` | 832 | 2.50e-3 | Fuel |
| `olive_oil` | 915 | 84.0e-3 | Viscous oil |
| `glycerin` | 1260 | 1.41 | Densest, by far the most viscous |

For anything else, pass `--density` (and optionally `--viscosity`), or set the
corresponding `MaterialConfig` fields, with the values for your liquid and
temperature.

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

- **Density drives the physics; liquid viscosity is carried but not yet coupled.**
  Viscosity is a first-class material property now (configurable per liquid,
  reported, and saved), but it does not yet change any trajectory. Note the
  Reynolds-dependent drag model uses the **air** viscosity `μ_air`
  (`Re = ρ_air·|v|·2r / μ_air`), *not* the liquid's — so the sprayed liquid's
  viscosity still has no effect on flight. In a fuller model it would govern
  **atomisation** — how the sheet breaks up and what droplet sizes form — and
  could feed the orifice discharge coefficient at low Reynolds numbers. Here the
  droplet size distribution is still an **input** (see the
  [sprayer parameters](sprayer_parameters.md)), so changing a liquid's viscosity
  alone leaves the result unchanged; changing its density does not. Liquid
  viscosity is in place as the hook for those future effects.
- **Surface tension** is not modelled at all.
- Properties are treated as constant (incompressible, isothermal). No temperature
  or evaporation effects.
- The material name, density, viscosity and solids fraction are saved into the
  `.npz` output, so every run is fully reproducible and self-describing. Archives
  written before a field existed still load, defaulting viscosity to water's and
  solids fraction to `1.0`.

See [sprayer_parameters.md](sprayer_parameters.md) for the other half of the
inputs.
