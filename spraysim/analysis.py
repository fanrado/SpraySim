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


@dataclass
class UniformityStats:
    """How even a deposited film is over a region of interest (ROI)."""

    n_cells: int
    mean_thickness: float        # m, over the ROI
    cv: float                    # coefficient of variation std/mean (0 = perfect)
    christiansen_cu: float       # Christiansen uniformity coefficient (1 = perfect)
    coverage_fraction: float     # fraction of ROI cells at/above the threshold
    coverage_threshold: float    # m, the threshold actually used
    min_thickness: float
    max_thickness: float
    p10_thickness: float
    p90_thickness: float

    def as_dict(self) -> dict:
        return {
            "n_cells": self.n_cells,
            "mean_thickness_um": self.mean_thickness * 1e6,
            "cv": self.cv,
            "christiansen_cu": self.christiansen_cu,
            "coverage_fraction": self.coverage_fraction,
            "coverage_threshold_um": self.coverage_threshold * 1e6,
            "min_thickness_um": self.min_thickness * 1e6,
            "max_thickness_um": self.max_thickness * 1e6,
            "p10_thickness_um": self.p10_thickness * 1e6,
            "p90_thickness_um": self.p90_thickness * 1e6,
        }


def _roi_mask(field: DepositionField,
              roi: tuple[float, float, float, float] | None) -> np.ndarray:
    """Boolean (ny, nx) mask of cells to include in a uniformity calculation.

    ``None`` -> the wetted (non-zero) cells; a rectangle ``(xmin, xmax, ymin,
    ymax)`` -> all cells (incl. dry ones) whose centres fall inside it, so gaps
    inside a target area are penalised.
    """
    if roi is None:
        return field.nonzero_mask()
    xmin, xmax, ymin, ymax = roi
    xc = 0.5 * (field.x_edges[:-1] + field.x_edges[1:])
    yc = 0.5 * (field.y_edges[:-1] + field.y_edges[1:])
    in_x = (xc >= xmin) & (xc <= xmax)
    in_y = (yc >= ymin) & (yc <= ymax)
    return np.outer(in_y, in_x)


def uniformity(
    field: DepositionField,
    *,
    roi: tuple[float, float, float, float] | None = None,
    coverage_threshold: float | None = None,
) -> UniformityStats:
    """Uniformity metrics of a deposited film over a region of interest.

    ``roi`` defaults to the wetted (non-zero) cells; pass a rectangle
    ``(xmin, xmax, ymin, ymax)`` to score a target area (dry gaps included).
    ``coverage_threshold`` (m) defaults to half the ROI mean thickness.

    Metrics: CV (std/mean, 0 = perfect), the Christiansen uniformity coefficient
    ``CU = 1 - Σ|xᵢ-x̄| / (n·x̄)`` (1 = perfect), and the coverage fraction (cells
    at or above the threshold), plus mean / min / max / p10 / p90.
    """
    values = field.thickness[_roi_mask(field, roi)]
    n = int(values.size)
    thr = 0.0 if coverage_threshold is None else float(coverage_threshold)

    if n == 0 or values.sum() <= 0.0:
        return UniformityStats(n, 0.0, 0.0, 0.0, 0.0, thr, 0.0, 0.0, 0.0, 0.0)

    mean = float(values.mean())
    cv = float(values.std()) / mean
    cu = 1.0 - float(np.abs(values - mean).sum()) / (n * mean)
    if coverage_threshold is None:
        thr = 0.5 * mean
    coverage = float(np.mean(values >= thr))
    return UniformityStats(
        n_cells=n,
        mean_thickness=mean,
        cv=cv,
        christiansen_cu=cu,
        coverage_fraction=coverage,
        coverage_threshold=thr,
        min_thickness=float(values.min()),
        max_thickness=float(values.max()),
        p10_thickness=float(np.percentile(values, 10)),
        p90_thickness=float(np.percentile(values, 90)),
    )


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
