"""Persist a :class:`SimResult` (and its config) to a ``.npz`` archive.

``numpy``'s ``.npz`` format stores a dictionary of named arrays, which is a
natural fit for the per-droplet result arrays. Two things need special care:

* **Trajectories** are a list of variable-length ``(steps, 3)`` arrays. Rather
  than an object array (which needs ``allow_pickle``), they are flattened into
  one ``(total_steps, 3)`` point array plus a per-trajectory ``lengths`` index,
  and split back apart on load.
* **Config** values are flattened into ``cfg_*`` scalars so the archive is
  self-describing and the exact run can be reconstructed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import PhysicsConfig, NozzleConfig, SimConfig
from .simulator import SimResult

# Bumped when the on-disk layout changes in an incompatible way.
FORMAT_VERSION = 1


def _flatten_config(cfg: SimConfig) -> dict[str, np.ndarray]:
    """Serialise a SimConfig into flat ``cfg_*`` entries (NaN encodes None)."""
    phys, noz = cfg.physics, cfg.nozzle
    n_dpl = np.nan if cfg.n_droplets is None else float(cfg.n_droplets)
    seed = np.nan if cfg.seed is None else float(cfg.seed)
    return {
        # Physics
        "cfg_gravity": phys.gravity,
        "cfg_air_density": phys.air_density,
        "cfg_water_density": phys.water_density,
        "cfg_drag_coefficient": phys.drag_coefficient,
        "cfg_ground_z": phys.ground_z,
        # Nozzle
        "cfg_position": np.asarray(noz.position, dtype=float),
        "cfg_direction": np.asarray(noz.direction, dtype=float),
        "cfg_half_angle": noz.half_angle,
        "cfg_pressure": noz.pressure,
        "cfg_orifice_diameter": noz.orifice_diameter,
        "cfg_shape": noz.shape,
        "cfg_speed_spread": noz.speed_spread,
        "cfg_distribution": noz.distribution,
        "cfg_mean_radius": noz.mean_radius,
        "cfg_radius_std": noz.radius_std,
        # Sim
        "cfg_n_droplets": n_dpl,
        "cfg_spray_duration": cfg.spray_duration,
        "cfg_max_droplets": cfg.max_droplets,
        "cfg_dt": cfg.dt,
        "cfg_max_time": cfg.max_time,
        "cfg_n_trajectories": cfg.n_trajectories,
        "cfg_seed": seed,
    }


def _rebuild_config(z: np.lib.npyio.NpzFile) -> SimConfig:
    """Inverse of :func:`_flatten_config`."""
    def f(key: str) -> float:
        return float(z[key])

    n_dpl = z["cfg_n_droplets"]
    seed = z["cfg_seed"]
    return SimConfig(
        n_droplets=None if np.isnan(n_dpl) else int(round(float(n_dpl))),
        spray_duration=f("cfg_spray_duration"),
        max_droplets=int(z["cfg_max_droplets"]),
        dt=f("cfg_dt"),
        max_time=f("cfg_max_time"),
        n_trajectories=int(z["cfg_n_trajectories"]),
        seed=None if np.isnan(seed) else int(round(float(seed))),
        physics=PhysicsConfig(
            gravity=f("cfg_gravity"),
            air_density=f("cfg_air_density"),
            water_density=f("cfg_water_density"),
            drag_coefficient=f("cfg_drag_coefficient"),
            ground_z=f("cfg_ground_z"),
        ),
        nozzle=NozzleConfig(
            position=tuple(z["cfg_position"].tolist()),
            direction=tuple(z["cfg_direction"].tolist()),
            half_angle=f("cfg_half_angle"),
            pressure=f("cfg_pressure"),
            orifice_diameter=f("cfg_orifice_diameter"),
            shape=str(z["cfg_shape"]),
            speed_spread=f("cfg_speed_spread"),
            distribution=str(z["cfg_distribution"]),
            mean_radius=f("cfg_mean_radius"),
            radius_std=f("cfg_radius_std"),
        ),
    )


def save_result(
    result: SimResult,
    config: SimConfig,
    path: str | Path,
    *,
    compress: bool = True,
) -> Path:
    """Write ``result`` (and ``config``) to a ``.npz`` archive at ``path``.

    The parent directory is created if needed. Set ``compress=False`` for a
    faster, larger uncompressed archive.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Flatten variable-length trajectories into points + lengths.
    if result.trajectories:
        traj_points = np.concatenate(result.trajectories, axis=0)
        traj_lengths = np.array([len(t) for t in result.trajectories], dtype=np.int64)
    else:
        traj_points = np.empty((0, 3), dtype=float)
        traj_lengths = np.empty(0, dtype=np.int64)

    data: dict[str, np.ndarray] = {
        "format_version": FORMAT_VERSION,
        # Per-droplet arrays.
        "landing_positions": result.landing_positions,
        "flight_times": result.flight_times,
        "impact_speeds": result.impact_speeds,
        "radii": result.radii,
        "launch_speeds": result.launch_speeds,
        "landed": result.landed,
        # Trajectories (flattened).
        "traj_points": traj_points,
        "traj_lengths": traj_lengths,
        # Derived hydraulics.
        "exit_speed": result.exit_speed,
        "flow_rate": result.flow_rate,
        "droplets_capped": result.droplets_capped,
    }
    data.update(_flatten_config(config))

    saver = np.savez_compressed if compress else np.savez
    saver(path, **data)
    # numpy appends .npz if missing; return the real path.
    if path.suffix != ".npz":
        path = path.with_name(path.name + ".npz")
    return path


def load_result(path: str | Path) -> tuple[SimResult, SimConfig]:
    """Load a ``.npz`` archive written by :func:`save_result`.

    Returns the reconstructed ``(SimResult, SimConfig)``.
    """
    path = Path(path)
    with np.load(path) as z:
        # Split the flattened trajectories back into a list.
        lengths = z["traj_lengths"]
        points = z["traj_points"]
        bounds = np.cumsum(lengths)[:-1]
        trajectories = [np.asarray(a) for a in np.split(points, bounds)] if lengths.size else []

        result = SimResult(
            landing_positions=z["landing_positions"],
            flight_times=z["flight_times"],
            impact_speeds=z["impact_speeds"],
            radii=z["radii"],
            launch_speeds=z["launch_speeds"],
            trajectories=trajectories,
            landed=z["landed"],
            exit_speed=float(z["exit_speed"]),
            flow_rate=float(z["flow_rate"]),
            droplets_capped=bool(z["droplets_capped"]),
        )
        config = _rebuild_config(z)
    return result, config
