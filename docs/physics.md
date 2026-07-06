# SpraySim — Physics Model

How a droplet is born, launched, and flown. This document collects the physics;
for the *inputs* that drive it see [sprayer_parameters.md](sprayer_parameters.md)
(the nozzle) and [material_properties.md](material_properties.md) (the liquid).

A nozzle emits droplets into a velocity cone; each droplet is a sphere of the
sprayed liquid, integrated under gravity and aerodynamic drag until it lands on
the ground plane.

## Forces and integration

Two forces act on each droplet:

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

Droplets do not interact, so a **moving nozzle** (a G-code toolpath) is just a
batch of droplets emitted from many positions along the path — the same
integrator handles it in one pass. A droplet also inherits the nozzle's
**carriage velocity**, so a fast-moving nozzle throws its droplets slightly
downstream (see [sprayer_parameters.md](sprayer_parameters.md)).

## Drag coefficient (Reynolds-dependent)

The drag coefficient `C_d` is not constant — it depends on the droplet
**Reynolds number**:

```
Re = ρ_air · |v| · (2r) / μ_air
```

The default `clift_gauvin` model recomputes `C_d(Re)` each step:

```
C_d(Re) = (24/Re)·(1 + 0.15·Re^0.687) + 0.42 / (1 + 4.25e4·Re^-1.16)
```

It reduces to **Stokes drag** (`C_d = 24/Re`) as `Re → 0` and to the **Newton**
plateau (~0.44) at high `Re`. This matters a lot for fine droplets: a fixed
`C_d = 0.47` (valid only for `Re ≳ 1000`) under-predicts drag at low `Re` and so
over-predicts terminal velocity, range and impact speed — by 2× or more below
~0.2 mm radius. A `constant` model (fixed `C_d`, the original behaviour) is
available via `--drag-model` / `DRAG_MODEL` for comparison and to reproduce older
runs.

Note the Reynolds number uses the **air** viscosity `μ_air`, not the liquid's;
the sprayed liquid's viscosity does not (yet) affect flight.

## Nozzle hydraulics (droplet count & exit speed)

Rather than typing a droplet count and speed directly, both are **derived** from
the nozzle pressure, orifice size and shape via incompressible orifice flow
(Bernoulli / Torricelli):

```
ideal velocity   v_ideal = √(2 ΔP / ρ_liquid)
exit speed       v       = C_v · v_ideal
volumetric flow  Q       = C_d · A · v_ideal          (A = orifice area)
droplet count    N       = Q · spray_duration / V_droplet
```

`C_d` (discharge) and `C_v` (velocity) coefficients come from the nozzle
**shape** (`sharp_orifice`, `rounded_orifice`, `full_cone`, `hollow_cone`,
`flat_fan`). `V_droplet = (4/3)π·E[r³]` is the mean droplet volume, so for a
fixed flow, finer droplets mean *many* more of them (`N ∝ 1/r³`).

> Note: `C_d` here is the **discharge** coefficient of the orifice — distinct
> from the aerodynamic **drag** coefficient `C_d(Re)` above.

## Droplet size distribution

Droplet radii are drawn from a configurable **`normal`** (Gaussian) or
**`lognormal`** distribution, both parameterised by a linear-space mean radius
and standard deviation. `E[r³]` (which sets the droplet count) is computed
analytically for each distribution:

- **normal:** the sampler clips radii at 0, so `E[r³]` is the partial
  (clipped-at-0) third moment, which matches the sampled droplets at any spread:
  `Φ(m/s)·(m³ + 3ms²) + φ(m/s)·s·(m² + 2s²)` (Φ, φ are the standard normal CDF and
  PDF). It reduces to `m³ + 3ms²` when `s ≪ m`. A wide normal (`m/s` small) also
  clips many radii to ~0, which is unphysical — the code warns and suggests
  `lognormal`.
- **lognormal:** `E[r³] = exp(3μ + 4.5σ²)`, with `σ² = ln(1 + (s/m)²)`,
  `μ = ln m − σ²/2`

## Material (sprayed liquid)

The simulation is not limited to water. The sprayed liquid's **density**
`ρ_liquid` enters the physics in two places:

- **Droplet mass** in the drag factor `k` above — a denser liquid gives heavier
  droplets that decelerate less (smaller `k`).
- **Exit hydraulics** — `v_ideal = √(2 ΔP / ρ_liquid)`, so a denser liquid exits
  more slowly for the same pressure (`v ∝ 1/√ρ_liquid`) and delivers a lower flow
  rate and droplet count.

These pull range in opposite directions (slower launch vs. less drag), which the
simulation resolves. Each material also carries a **dynamic viscosity**, but that
is not yet coupled to the flight physics — see
[material_properties.md](material_properties.md) for the material registry and
the full discussion.
