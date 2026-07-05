"""Derived statistics from a simulation result."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import SimConfig
from .simulator import SimResult


@dataclass
class SprayStats:
    n_droplets: int
    landed_fraction: float
    mean_flight_time: float
    median_flight_time: float
    mean_impact_speed: float
    mean_radius_mm: float
    coverage_radius_p50: float   # median landing distance from the spray axis
    coverage_radius_p90: float   # 90th-percentile landing distance
    spray_center: tuple[float, float]

    def as_dict(self) -> dict:
        return {
            "n_droplets": self.n_droplets,
            "landed_fraction": self.landed_fraction,
            "mean_flight_time_s": self.mean_flight_time,
            "median_flight_time_s": self.median_flight_time,
            "mean_impact_speed_ms": self.mean_impact_speed,
            "mean_radius_mm": self.mean_radius_mm,
            "coverage_radius_p50_m": self.coverage_radius_p50,
            "coverage_radius_p90_m": self.coverage_radius_p90,
            "spray_center_xy_m": self.spray_center,
        }


def radial_distances(result: SimResult, config: SimConfig) -> np.ndarray:
    """Horizontal distance of each landing point from the nozzle's x,y."""
    cx, cy, _ = config.nozzle.position
    dx = result.landing_positions[:, 0] - cx
    dy = result.landing_positions[:, 1] - cy
    return np.hypot(dx, dy)


@dataclass
class DepositionField:
    """A 2-D map of dry deposited film thickness on the ground plane.

    ``thickness`` is indexed ``[iy, ix]`` (rows = y, cols = x) to match image /
    ``imshow`` conventions, in metres of *dry* (cured) film.
    """

    x_edges: np.ndarray          # (nx+1,) cell boundaries in x (m)
    y_edges: np.ndarray          # (ny+1,) cell boundaries in y (m)
    thickness: np.ndarray        # (ny, nx) dry film thickness (m)
    cell_size: float             # m, square cell edge
    solids_fraction: float       # volume fraction of solids used for the dry film

    @property
    def cell_area(self) -> float:
        return self.cell_size ** 2

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """(xmin, xmax, ymin, ymax) — e.g. for ``imshow(extent=...)``."""
        return (float(self.x_edges[0]), float(self.x_edges[-1]),
                float(self.y_edges[0]), float(self.y_edges[-1]))

    def nonzero_mask(self) -> np.ndarray:
        """Boolean mask of cells that received any deposit (the wetted region)."""
        return self.thickness > 0.0

    def total_dry_volume(self) -> float:
        return float(self.thickness.sum()) * self.cell_area

    def total_wet_volume(self) -> float:
        sf = self.solids_fraction
        return self.total_dry_volume() / sf if sf > 0.0 else 0.0


def deposition_map(
    result: SimResult,
    config: SimConfig,
    *,
    cell_size: float | None = None,
    extent: tuple[float, float, float, float] | None = None,
) -> DepositionField:
    """Bin landed droplets into a dry film-thickness map.

    Each landed droplet deposits its volume ``(4/3)π r³`` where it lands; the wet
    thickness of a cell is deposited volume / cell area, and the dry thickness is
    that times ``config.material.solids_fraction`` (no spreading/run-off/
    evaporation kinetics are modelled — deposit-where-it-lands).

    ``cell_size`` (m) defaults to ~1/60 of the larger landing extent. ``extent``
    (xmin, xmax, ymin, ymax) defaults to the bounding box of the landings;
    droplets outside a supplied extent are ignored.
    """
    landed = result.landed
    pos = result.landing_positions[landed]
    radii = result.radii[landed]

    if extent is None:
        if pos.shape[0] == 0:
            xmin = xmax = ymin = ymax = 0.0
        else:
            xmin, ymin = float(pos[:, 0].min()), float(pos[:, 1].min())
            xmax, ymax = float(pos[:, 0].max()), float(pos[:, 1].max())
    else:
        xmin, xmax, ymin, ymax = (float(v) for v in extent)

    if cell_size is None:
        span = max(xmax - xmin, ymax - ymin)
        cell_size = span / 60.0 if span > 0.0 else 1.0e-3
    cell_size = float(cell_size)

    nx = max(1, int(np.ceil((xmax - xmin) / cell_size)))
    ny = max(1, int(np.ceil((ymax - ymin) / cell_size)))
    x_edges = xmin + cell_size * np.arange(nx + 1)
    y_edges = ymin + cell_size * np.arange(ny + 1)

    sf = config.material.solids_fraction
    thickness = np.zeros((ny, nx))
    if pos.shape[0] > 0:
        dry_volume = (4.0 / 3.0) * np.pi * radii ** 3 * sf
        # histogram2d returns (nx, ny) indexed [ix, iy]; transpose to [iy, ix].
        binned, _, _ = np.histogram2d(
            pos[:, 0], pos[:, 1], bins=[x_edges, y_edges], weights=dry_volume
        )
        thickness = binned.T / (cell_size ** 2)

    return DepositionField(x_edges, y_edges, thickness, cell_size, sf)


def summarize(result: SimResult, config: SimConfig) -> SprayStats:
    landed = result.landed
    radial = radial_distances(result, config)
    landed_radial = radial[landed] if landed.any() else radial

    return SprayStats(
        n_droplets=result.n,
        landed_fraction=float(np.mean(landed)),
        mean_flight_time=float(np.mean(result.flight_times)),
        median_flight_time=float(np.median(result.flight_times)),
        mean_impact_speed=float(np.mean(result.impact_speeds)),
        mean_radius_mm=float(np.mean(result.radii) * 1000.0),
        coverage_radius_p50=float(np.percentile(landed_radial, 50)),
        coverage_radius_p90=float(np.percentile(landed_radial, 90)),
        spray_center=(
            float(np.mean(result.landing_positions[landed, 0])) if landed.any() else 0.0,
            float(np.mean(result.landing_positions[landed, 1])) if landed.any() else 0.0,
        ),
    )
