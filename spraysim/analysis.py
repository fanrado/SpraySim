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
