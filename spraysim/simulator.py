"""Core integrator: advances all droplets until they land or time out."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import SimConfig
from .nozzle import Nozzle


@dataclass
class SimResult:
    """Outcome of a simulation run. All arrays are ordered by droplet index."""

    landing_positions: np.ndarray   # (n, 3) where each droplet hit the ground
    flight_times: np.ndarray        # (n,) seconds aloft
    impact_speeds: np.ndarray       # (n,) speed magnitude at impact (m/s)
    radii: np.ndarray               # (n,) droplet radius (m)
    launch_speeds: np.ndarray       # (n,) initial speed (m/s)
    trajectories: list[np.ndarray]  # sampled full paths, each (steps, 3)
    landed: np.ndarray              # (n,) bool, False if it timed out mid-air

    @property
    def n(self) -> int:
        return self.landing_positions.shape[0]


class Simulator:
    """Semi-implicit Euler integrator with quadratic aerodynamic drag."""

    def __init__(self, config: SimConfig | None = None):
        self.config = config or SimConfig()

    def run(self) -> SimResult:
        cfg = self.config
        phys = cfg.physics
        rng = np.random.default_rng(cfg.seed)

        nozzle = Nozzle(cfg.nozzle)
        pos, vel, radii = nozzle.emit(cfg.n_droplets, rng)
        launch_speeds = np.linalg.norm(vel, axis=1)

        # Per-droplet drag factor k so that a_drag = -k * |v| * v.
        # Derived from F_drag = 0.5 * rho_air * Cd * A * |v| * v with
        # A = pi r^2 and m = rho_water * (4/3) pi r^3.
        k = (3.0 * phys.air_density * phys.drag_coefficient) / (
            8.0 * phys.water_density * radii
        )
        gravity = np.array([0.0, 0.0, -phys.gravity])

        n = cfg.n_droplets
        active = np.ones(n, dtype=bool)
        time_aloft = np.zeros(n)

        landing_positions = pos.copy()
        impact_speeds = np.zeros(n)
        landed = np.zeros(n, dtype=bool)

        # Record full paths for a deterministic sample of droplets.
        sample_idx = np.linspace(
            0, n - 1, min(cfg.n_trajectories, n), dtype=int
        )
        sample_set = set(sample_idx.tolist())
        traj_store: dict[int, list] = {i: [pos[i].copy()] for i in sample_set}

        dt = cfg.dt
        max_steps = int(np.ceil(cfg.max_time / dt))

        for _ in range(max_steps):
            if not active.any():
                break

            a = active.copy()  # snapshot: active is mutated as droplets land
            speed = np.linalg.norm(vel[a], axis=1)
            accel = gravity - (k[a] * speed)[:, None] * vel[a]

            # Semi-implicit (symplectic) Euler: update velocity, then position.
            vel[a] += accel * dt
            new_pos = pos[a] + vel[a] * dt

            # Detect ground crossing and clamp to the exact impact point.
            below = new_pos[:, 2] <= phys.ground_z
            act_idx = np.nonzero(a)[0]

            if below.any():
                z0 = pos[a][below, 2]
                z1 = new_pos[below, 2]
                # Linear fraction of the step at which z hits ground_z.
                frac = (z0 - phys.ground_z) / np.clip(z0 - z1, 1e-12, None)
                frac = np.clip(frac, 0.0, 1.0)
                hit = pos[a][below] + frac[:, None] * (new_pos[below] - pos[a][below])
                hit[:, 2] = phys.ground_z

                landed_idx = act_idx[below]
                landing_positions[landed_idx] = hit
                impact_speeds[landed_idx] = np.linalg.norm(vel[a][below], axis=1)
                time_aloft[landed_idx] += frac * dt
                landed[landed_idx] = True
                active[landed_idx] = False

            pos[a] = new_pos
            time_aloft[act_idx[~below]] += dt

            # Append sampled trajectory points for still-active sampled droplets.
            for i in sample_set:
                if active[i]:
                    traj_store[i].append(pos[i].copy())

        # Any droplet still active timed out; record its last known state.
        still = np.nonzero(active)[0]
        if still.size:
            landing_positions[still] = pos[still]
            impact_speeds[still] = np.linalg.norm(vel[still], axis=1)

        trajectories = [np.asarray(traj_store[i]) for i in sample_idx.tolist()]

        return SimResult(
            landing_positions=landing_positions,
            flight_times=time_aloft,
            impact_speeds=impact_speeds,
            radii=radii,
            launch_speeds=launch_speeds,
            trajectories=trajectories,
            landed=landed,
        )
