# Sprayer Parameters

This is one of the **two inputs that dictate a spray result** (the other is the
[material](material_properties.md)). Everything here describes the *sprayer* —
where it is, how it is aimed, how hard it pushes the liquid, and what droplets it
produces. Together with the liquid's density these parameters fully determine the
droplet count, exit speed, spread and landing pattern.

All quantities are **SI units** (metres, seconds, kilograms, pascals, radians)
inside the code. The config files and CLI accept friendlier units (bar, mm,
degrees) and convert them.

- Config keys live in `config/*.conf` (sourced by `main.sh`).
- CLI flags are passed to `run.py`.
- Library fields live on `NozzleConfig` / `SimConfig` in `spraysim/config.py`.

---

## 1. Geometry — where droplets are born and where they go

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Height | `HEIGHT` | `--height` | `NozzleConfig.position[2]` | m | `1.5` | Nozzle height above the ground plane. Higher → longer flight, wider landing pattern. |
| Position | — | — | `NozzleConfig.position` | m | `(0, 0, 1.5)` | Full 3-D nozzle location. Only settable via the library; the CLI exposes height. |
| Direction | — | — | `NozzleConfig.direction` | unit vector | `(0, 0, -1)` | Spray axis (need not be normalised). Default points straight down. |
| Cone half-angle | `CONE` | `--cone` | `NozzleConfig.half_angle` | degrees (CLI) / radians (field) | `25°` | Half of the full cone spread. Droplets are emitted uniformly over this solid angle. Wider cone → larger, more diffuse footprint. |

**How the cone is sampled:** directions are drawn with `cos θ` uniform in
`[cos(half_angle), 1]`, which spreads droplets *evenly per unit area* on the cone
cap (no bunching at the axis).

---

## 2. Hydraulics — pressure and orifice set the count and speed

Rather than typing a droplet count and speed, both are **derived** from the
pressure, orifice size, and nozzle shape via incompressible orifice flow
(Bernoulli / Torricelli). See `spraysim/hydraulics.py`.

```
ideal velocity   v_ideal = √(2 · ΔP / ρ_liquid)
exit speed       v       = C_v · v_ideal
volumetric flow  Q       = C_d · A · v_ideal          (A = π (d/2)²)
droplet count    N       = Q · spray_duration / V_droplet
```

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Pressure | `PRESSURE_BAR` | `--pressure-bar` | `NozzleConfig.pressure` | bar (CLI) / Pa (field) | `3.0` bar | Gauge pressure drop across the orifice. Drives exit speed and flow. |
| Orifice diameter | `ORIFICE_MM` | `--orifice-mm` | `NozzleConfig.orifice_diameter` | mm (CLI) / m (field) | `0.8` mm | Hole diameter. Flow ∝ area ∝ diameter², so a big lever on droplet count. |
| Nozzle shape | `NOZZLE_SHAPE` | `--shape` | `NozzleConfig.shape` | name | `full_cone` | Selects the discharge/velocity coefficients (table below). |
| Spray duration | `SPRAY_DURATION` | `--spray-duration` | `SimConfig.spray_duration` | s | `0.15` | How long the nozzle is open. Linearly scales the droplet count. |
| Explicit count | `DROPLETS` | `--droplets` | `SimConfig.n_droplets` | integer | *(empty → derive)* | Optional override that pins the droplet count and skips the derivation. |
| Count cap | — | — | `SimConfig.max_droplets` | integer | `200000` | Safety cap on the derived count; a run reports if it was hit. |

### Nozzle shapes → coefficients

`C_d` (discharge) folds in the vena-contracta contraction; `C_v` (velocity) is
the velocity loss. Higher `C_d` → more flow (more droplets); higher `C_v` →
faster droplets.

| Shape | `C_d` | `C_v` | Typical use |
|-------|-------|-------|-------------|
| `sharp_orifice` | 0.61 | 0.98 | Thin, sharp-edged hole |
| `rounded_orifice` | 0.92 | 0.98 | Well-rounded inlet, little contraction |
| `full_cone` | 0.75 | 0.95 | Solid-cone agricultural nozzle *(default)* |
| `hollow_cone` | 0.70 | 0.95 | Swirl / hollow-cone atomiser |
| `flat_fan` | 0.88 | 0.97 | Flat-fan spray tip |

### What moves the result

| Increase… | Exit speed `v` | Flow `Q` | Droplet count `N` |
|-----------|----------------|----------|-------------------|
| Pressure `ΔP` | ↑ (∝ √ΔP) | ↑ (∝ √ΔP) | ↑ (∝ √ΔP) |
| Orifice diameter `d` | — | ↑ (∝ d²) | ↑ (∝ d²) |
| Spray duration | — | — | ↑ (∝ duration) |
| Droplet size (mean radius) | — | — | ↓ (∝ 1/r³) — see below |

---

## 3. Droplet size distribution

Droplet radii are drawn from a configurable distribution. The **mean cubed
radius** `E[r³]` sets the mean droplet volume `V_droplet = (4/3)π·E[r³]`, which in
turn sets the droplet count (`N = Q·duration / V_droplet`). Because volume goes
as `r³`, finer droplets mean *many* more of them for the same flow.

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Distribution | `DISTRIBUTION` | `--distribution` | `NozzleConfig.distribution` | `normal` \| `lognormal` | `lognormal` | Shape of the radius distribution. |
| Mean radius | `MEAN_RADIUS_MM` | `--mean-radius-mm` | `NozzleConfig.mean_radius` | mm (CLI) / m (field) | `0.4` mm | Linear-space mean droplet radius. |
| Radius std | `RADIUS_STD_MM` | `--radius-std-mm` | `NozzleConfig.radius_std` | mm (CLI) / m (field) | `0.12` mm | Linear-space standard deviation of the radius. |

Both distributions are parameterised by the **linear-space mean and std**; the
model converts these to the underlying parameters and uses a closed form for
`E[r³]`:

- **normal:** the clipped-at-0 partial third moment
  `Φ(m/s)·(m³+3ms²) + φ(m/s)·s·(m²+2s²)` (matches the sampled droplets at any
  spread; reduces to `m³+3ms²` for `s ≪ m`)
- **lognormal:** `E[r³] = exp(3μ + 4.5σ²)`, with `σ² = ln(1 + (s/m)²)`, `μ = ln m − σ²/2`

Use `normal` for a tight, symmetric mist and `lognormal` for a realistic
right-skewed spray with a tail of larger drops. A wide `normal` (`radius_std`
comparable to `mean_radius`) clips many radii to ~0 and **warns** — prefer
`lognormal` there.

---

## 4. Emission speed spread

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Speed spread | `SPEED_SPREAD` | `--speed-spread` | `NozzleConfig.speed_spread` | fraction | `0.15` | Relative std of per-droplet speed about the derived exit speed (turbulence). `0` = every droplet leaves at exactly the exit speed. |

Speeds are drawn from `Normal(exit_speed, speed_spread · exit_speed)` and clipped
to be non-negative.

---

## 5. Simulation / integration controls

These do not change the physics; they control accuracy, runtime, and
reproducibility.

| Parameter | Config key | CLI flag | Field | Unit | Default | Meaning |
|-----------|-----------|----------|-------|------|---------|---------|
| Timestep | `DT` | `--dt` | `SimConfig.dt` | s | `0.001` | Integration step. Smaller → more accurate, slower. |
| Max flight time | — | — | `SimConfig.max_time` | s | `8.0` | Hard cap; any droplet still aloft is recorded at its last state. |
| Sampled trajectories | — | — | `SimConfig.n_trajectories` | count | `60` | Number of droplets whose full path is stored for the side-view plot. |
| Seed | `SEED` | `--seed` | `SimConfig.seed` | integer | `42` | RNG seed for reproducible runs. |

---

## 6. Environment constants (not the sprayer, but they act on the spray)

Held on `PhysicsConfig`. `drag_model` is exposed as a config key / CLI flag; the
rest are library-only but influence every trajectory.

| Field | Config key | CLI flag | Unit | Default | Meaning |
|-------|-----------|----------|------|---------|---------|
| `gravity` | — | — | m/s² | `9.81` | Downward acceleration. |
| `air_density` | — | — | kg/m³ | `1.225` | Sets the drag magnitude (sea level, 15 °C). `0` disables drag. |
| `air_viscosity` | — | — | Pa·s | `1.81e-5` | Air dynamic viscosity; sets the droplet Reynolds number. |
| `drag_model` | `DRAG_MODEL` | `--drag-model` | name | `clift_gauvin` | `clift_gauvin` (Reynolds-dependent `C_d`) or `constant` (fixed `C_d`, legacy). |
| `drag_coefficient` | — | — | — | `0.47` | Fixed `C_d` used **only** by the `constant` drag model. |
| `ground_z` | — | — | m | `0.0` | Height of the impact plane. |

The per-droplet drag factor is `k = 3·ρ_air·C_d / (8·ρ_liquid·r)`, so the
**air density**, the **liquid density from the material**, and — through
`C_d(Re)` with `Re = ρ_air·|v|·2r / μ_air` — the **air viscosity** all set how
quickly droplets slow down.

The `clift_gauvin` model recomputes `C_d(Re)` every step; it reduces to Stokes
drag at low `Re` and the Newton plateau at high `Re`. A fixed `C_d` (the
`constant` model) over-predicts range and impact speed for fine droplets by 2× or
more — see [material_properties.md](material_properties.md).
