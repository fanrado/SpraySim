"""Core integrator: advances all droplets until they land or time out."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import SimConfig
from .nozzle import Nozzle
from . import hydraulics, drag, gcode


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
    # Derived nozzle hydraulics (constant for the run).
    exit_speed: float = 0.0         # m/s, pressure-derived launch speed
    flow_rate: float = 0.0          # m^3/s, volumetric flow through the orifice
    droplets_capped: bool = False   # True if the derived count hit max_droplets
    # Spray-on (G1) segments of the toolpath, (k, 4) [x0, y0, x1, y1]; None if
    # the run sprayed from a single fixed position.
    path_segments: np.ndarray | None = None

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

        nozzle = Nozzle(cfg.nozzle, cfg.material.density)

        # Droplet count: explicit override, or derived from flow * spray time.
        # For a G-code path the "spray time" is the total time on G1 moves.
        droplets_capped = False
        path_segments = None
        if cfg.path is not None:
            moves = gcode.load_moves(
                cfg.path.gcode, default_feed=cfg.path.default_feed,
                standoff=cfg.path.standoff, feed_override=cfg.path.feed_override,
            )
            spray_time = gcode.total_spray_time(moves)
            path_segments = np.array(
                [[m.start[0], m.start[1], m.end[0], m.end[1]]
                 for m in moves if m.spray_on], dtype=float
            )
        else:
            moves = None
            spray_time = cfg.spray_duration

        if cfg.n_droplets is not None:
            n = int(cfg.n_droplets)
        else:
            raw_n = hydraulics.droplet_count(
                nozzle.flow_rate, spray_time, nozzle.mean_droplet_volume()
            )
            n = min(raw_n, cfg.max_droplets)
            droplets_capped = raw_n > cfg.max_droplets

        if moves is not None:
            pos, vel, radii = nozzle.emit_path(
                moves, n, rng,
                include_carriage_velocity=cfg.path.include_carriage_velocity,
            )
            n = pos.shape[0]
            if n == 0:
                raise ValueError("G-code path has no spray (G1) moves to emit from.")
        else:
            pos, vel, radii = nozzle.emit(n, rng)
        launch_speeds = np.linalg.norm(vel, axis=1)

        # Per-droplet drag factor k so that a_drag = -k * |v| * v, derived from
        # F_drag = 0.5 * rho_air * Cd * A * |v| * v with A = pi r^2 and
        # m = rho_liquid * (4/3) pi r^3, i.e. k = 3 rho_air Cd / (8 rho_liquid r).
        #
        # "constant" model: Cd is fixed, so k is precomputed once (identical to
        # the original arithmetic). "clift_gauvin": Cd = Cd(Re) is recomputed per
        # step, so only the radius-dependent base factor is precomputed here.
        constant_drag = phys.drag_model == drag.CONSTANT
        if constant_drag:
            k_const = (3.0 * phys.air_density * phys.drag_coefficient) / (
                8.0 * cfg.material.density * radii
            )
        else:
            k_base = (3.0 * phys.air_density) / (8.0 * cfg.material.density * radii)
        gravity = np.array([0.0, 0.0, -phys.gravity])

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
            pos_a = pos[a]
            v_before = vel[a]  # velocity at the start of the step (copy)
            speed = np.linalg.norm(v_before, axis=1)

            # Drag factor for the active droplets this step.
            if constant_drag:
                k_a = k_const[a]
            else:
                Re = drag.reynolds_number(speed, radii[a], phys.air_density,
                                          phys.air_viscosity)
                k_a = k_base[a] * drag.drag_coefficient(Re, phys.drag_model)

            accel = gravity - (k_a * speed)[:, None] * v_before

            # Semi-implicit (symplectic) Euler: update velocity, then position.
            v_after = v_before + accel * dt
            new_pos = pos_a + v_after * dt

            # Detect ground crossing and clamp to the exact impact point.
            below = new_pos[:, 2] <= phys.ground_z
            act_idx = np.nonzero(a)[0]

            if below.any():
                z0 = pos_a[below, 2]
                z1 = new_pos[below, 2]
                # Linear fraction of the step at which z hits ground_z.
                frac = (z0 - phys.ground_z) / np.clip(z0 - z1, 1e-12, None)
                frac = np.clip(frac, 0.0, 1.0)
                hit = pos_a[below] + frac[:, None] * (new_pos[below] - pos_a[below])
                hit[:, 2] = phys.ground_z
                # Impact velocity at the crossing: interpolate between the pre-
                # and post-step velocities at the same fraction as the position.
                # (Reading v_after would sample a partial step too late.)
                v_hit = v_before[below] + frac[:, None] * (v_after[below] - v_before[below])

                landed_idx = act_idx[below]
                landing_positions[landed_idx] = hit
                impact_speeds[landed_idx] = np.linalg.norm(v_hit, axis=1)
                time_aloft[landed_idx] += frac * dt
                landed[landed_idx] = True
                active[landed_idx] = False

            vel[a] = v_after
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
            exit_speed=nozzle.exit_speed,
            flow_rate=nozzle.flow_rate,
            droplets_capped=droplets_capped,
            path_segments=path_segments,
        )
