# SpraySim — Physics Fix Plan

Planned fixes from the physics validation pass. The integrator, drag law,
hydraulics and size-distribution moments were all verified correct to
discretization precision (terminal velocity to 1e-13, clean O(dt) convergence).
The items below address **model fidelity** and two minor numerical details, not
implementation bugs.

Priority: **P1** materially changes results and is the reason fine-mist runs are
unphysical; **P2/P3** are small-magnitude corrections; **T** is supporting
tooling.

---

## P1 — Reynolds-dependent drag coefficient  ✅ DONE

Implemented: `spraysim/drag.py` (Clift–Gauvin + constant models), `air_viscosity`
and `drag_model` on `PhysicsConfig` (default `clift_gauvin`), per-step `Cd(Re)` in
the integrator (constant path kept bit-for-bit), storage round-trip with old
archives falling back to `constant`, `--drag-model` / `DRAG_MODEL` plumbing, tests
and docs. Verified: constant reproduces legacy `v_t`; clift_gauvin matches the
Cd(Re) fixed-point solve; vacuum stays drag-free.

### Problem
`spraysim/simulator.py:65` builds the drag factor with a **constant**
`Cd = 0.47` (`PhysicsConfig.drag_coefficient`). That value is only valid in the
Newton regime (Re ≳ 1000). Real droplets frequently sit well below it, where drag
is far stronger, so the sim **over-predicts terminal velocity, range and impact
speed for fine droplets**.

Measured overestimate of terminal velocity vs a Reynolds-dependent correlation:

| radius | Re @ v_t | overestimate |
|--------|----------|--------------|
| 0.10 mm | 29 | +201% |
| 0.40 mm | 232 | +35% |
| 0.80 mm | 657 | +5% |
| 1.60 mm | 1857 | −3% |

Accurate for ≥0.8 mm; badly off for the whole `fine_mist` regime.

### Root cause
Drag coefficient does not depend on the droplet Reynolds number
`Re = ρ_air · |v| · d / μ_air` (d = 2r). Air viscosity `μ_air` is not modelled at
all.

### Proposed fix
Compute `Cd = Cd(Re)` per droplet, per step, using a smooth correlation valid
across the full range — **Clift–Gauvin**:

```
Cd(Re) = (24/Re)·(1 + 0.15·Re^0.687) + 0.42 / (1 + 4.25e4·Re^-1.16)
```

- Reduces to **Stokes** (`Cd = 24/Re`) as Re→0, which makes the drag term
  linear in v and matches the analytic Stokes acceleration `9μ_air v / (2 ρ_l r²)`
  (verified by hand — the limit is exact).
- Approaches ~0.44 (Newton) at high Re, so large-droplet behaviour is unchanged.
- Guard `v = 0` / `Re = 0` → zero drag (no division blowup).

Per step the acceleration becomes:
```
Re     = ρ_air · |v| · 2r / μ_air
Cd     = Cd(Re)
a_drag = -(3 ρ_air Cd / (8 ρ_l r)) · |v| · v
```
i.e. the existing `k` gains a per-step `Cd(Re)` factor instead of a constant.

### Design decisions
- Add `air_viscosity` to `PhysicsConfig` (default `1.81e-5` Pa·s, air at ~15 °C).
- Add a `drag_model` selector: `"clift_gauvin"` (new default) and `"constant"`
  (uses the existing fixed `drag_coefficient`, retained for A/B comparison and to
  reproduce old runs). Put the `Cd(Re)` correlations in a new module
  `spraysim/drag.py` (parallels `hydraulics.py` / `materials.py`).
- Keep `drag_coefficient` on `PhysicsConfig` — it is the constant used by the
  `"constant"` model.

### Files
- **new** `spraysim/drag.py` — `reynolds_number(...)`, `drag_coefficient(Re, model)`, model registry.
- `spraysim/config.py` — `PhysicsConfig`: add `air_viscosity`, `drag_model`.
- `spraysim/simulator.py` — recompute `Cd(Re)` inside the step loop (vectorised over active droplets); build `k` from it.
- `spraysim/storage.py` — persist/restore `cfg_air_viscosity`, `cfg_drag_model` (with fallback defaults for older archives, as done for viscosity).
- `spraysim/__init__.py` — export `drag`.
- `run.py` / `main.sh` / `config/*.conf` — optional `--drag-model` flag + `DRAG_MODEL` key (default clift_gauvin).
- `docs/sprayer_parameters.md` (environment section) + `README.md` — document the drag model and air viscosity.

### Acceptance criteria
- New test: simulated terminal velocity matches a fixed-point solve of
  `g = k(Cd(Re))·v²` to <1% across r ∈ [0.05, 1.6] mm.
- Stokes limit: for a very small droplet, terminal velocity matches
  `v_t = 2 r² (ρ_l−ρ_air) g / (9 μ_air)` (Stokes) to <2%.
- `drag_model="constant"` reproduces current results bit-for-bit.
- Existing tests still pass (vacuum unaffected; landing/stats qualitative).

### Risks
- Per-step `Cd(Re)` adds cost; keep it vectorised over active droplets only.
- Changes default results — call this out in the changelog/README; `"constant"`
  preserves the old behaviour.

---

## P2 — Normal-distribution `E[r³]` clipping bias

### Problem
`spraysim/nozzle.py:74` returns the analytic third moment `m³ + 3ms²` for the
`normal` distribution, but `sample_radii` clips `r ≤ 0`. For wide distributions
the analytic value drifts from the sampled droplets, so the derived droplet count
is off:

| s/m | E[r³] error |
|-----|-------------|
| 0.25 (default) | −0.06% |
| 1.0 | −1.9% |
| 1.5 | −7.9% |

Negligible for shipped presets (s/m ≈ 0.25–0.33); a latent trap for wide inputs.

### Proposed fix (pick one)
1. **Guardrail (low effort):** warn (or raise) when `s/m` exceeds a threshold
   (e.g. 0.4) for the normal distribution, steering users to `lognormal`.
2. **Truncated-moment correction (medium):** compute `E[r³]` for the
   left-truncated (at 0) normal in closed form so it matches the clipped sampler.
3. **Nudge default presets** toward `lognormal` where skew is expected.

Recommend **(1)** now (cheap, honest) and note (2) as a follow-up. Keep the
existing "valid while std ≪ mean" comment.

### Files
- `spraysim/nozzle.py` (validation/warn or truncated moment).
- `tests/test_simulator.py` — extend the moment test to a wide-`s/m` case.

### Acceptance criteria
- Either a clear warning fires for `s/m > threshold`, or analytic vs sampled
  `E[r³]` agree to <1% up to `s/m = 1`.

---

## P3 — Impact speed evaluated at end-of-step, not at ground crossing  ✅ DONE

Implemented in `spraysim/simulator.py`: snapshot the pre-step velocity and
interpolate the impact velocity to the same `frac` used for the landing position,
so speed and position are reported at the *same* crossing instant (previously
velocity was read a full step later). This is primarily a **consistency** fix —
whether it reduces error vs a fixed analytic value is scenario-dependent because
semi-implicit Euler has its own O(dt) position bias, but it removes the one-step
lag and, in drag-dominated cases, tightens impact speed vs the dt→0 reference
(~0.30% → ~0.23% at dt=1e-3). Validated by a vacuum energy-conservation test:
impact speed matches `sqrt(u²+2gh)` to 0.0033% at dt=1e-3, converging O(dt).
Non-landing trajectories are unchanged (constant-model reproduction preserved).


### Problem
`spraysim/simulator.py:114` interpolates the landing **position** to the exact
ground crossing (good) but reads the **velocity** at the end of the step. This is
an O(dt) inconsistency:

| dt | impact speed | landing x |
|----|--------------|-----------|
| 2e-3 | 6.699 | 5.850 |
| 1e-3 | 6.711 | 5.873 |
| 2e-6 (ref) | 6.731 | 5.896 |

~0.3% low on speed, ~0.4% on range at the default `dt = 1e-3`.

### Proposed fix
Evaluate impact speed at the interpolated crossing: use the pre-step velocity and
the same `frac` used for the position, i.e. blend `vel_before` and `vel_after`
(or recompute acceleration once at the crossing). Also advance `time_aloft` and
the landing position consistently with that `frac` (position already is).

### Files
- `spraysim/simulator.py` — snapshot `vel` before the update; interpolate at `frac`.
- `tests/test_simulator.py` — assert impact speed converges and the dt-sensitivity shrinks.

### Acceptance criteria
- Impact speed at `dt = 1e-3` within ~0.05% of the `dt → 0` reference (down from ~0.3%).

---

## T — Lock in the validation harness

The checks above were run from throwaway scripts. Make validation repeatable.

### Proposed
- Add `analysis/validate.py`: terminal velocity, vacuum free-fall, dt convergence,
  drag monotonicity, hydraulics identities, and (after P1) the Re-drag / Stokes
  benchmarks — printing a pass/fail table.
- Promote the cheap, deterministic ones into `tests/` as assertions
  (terminal velocity, convergence order, Stokes limit).

### Files
- **new** `analysis/validate.py`; extend `tests/test_simulator.py`;
  mention in `analysis/README.md`.

---

## Suggested sequencing

1. **T (harness first)** — so P1's before/after is measurable and guarded.
2. **P1 (Reynolds drag)** — the substantive fidelity fix; also gives `μ_air`
   (and, conceptually, viscosity) a real role.
3. **P3** then **P2** — small corrections, quick once P1's tests exist.

Each step: implement → run `pytest` + `analysis/validate.py` → commit file-by-file
with descriptive messages (per the project's commit convention).
