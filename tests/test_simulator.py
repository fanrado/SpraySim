"""Sanity tests for the spray simulator's physics and bookkeeping."""

import math

import numpy as np
import pytest

from spraysim import SimConfig, NozzleConfig, PhysicsConfig, Simulator, analysis
from spraysim.nozzle import Nozzle


def test_all_droplets_land_within_bounds():
    cfg = SimConfig(n_droplets=500, seed=1)
    result = Simulator(cfg).run()

    assert result.n == 500
    # With a downward spray and 8 s budget, every droplet should land.
    assert result.landed.all()
    # Nothing ends below the ground plane.
    assert np.all(result.landing_positions[:, 2] >= cfg.physics.ground_z - 1e-6)
    # Flight times are positive and finite.
    assert np.all(result.flight_times > 0)
    assert np.all(np.isfinite(result.flight_times))


def test_vacuum_freefall_matches_analytic():
    """With no drag and no launch velocity, z(t) = h - 0.5 g t^2."""
    h = 5.0
    phys = PhysicsConfig(air_density=0.0)  # disable drag
    nozzle = NozzleConfig(position=(0.0, 0.0, h), half_angle=0.0, exit_speed=0.0,
                          speed_spread=0.0)
    cfg = SimConfig(n_droplets=1, dt=1e-4, nozzle=nozzle, physics=phys, seed=0)

    result = Simulator(cfg).run()
    expected = math.sqrt(2 * h / phys.gravity)
    assert result.flight_times[0] == pytest.approx(expected, rel=1e-2)


def test_cone_directions_respect_half_angle():
    half = math.radians(20.0)
    cfg = NozzleConfig(direction=(0, 0, -1), half_angle=half, speed_spread=0.0)
    rng = np.random.default_rng(0)
    _, vel, _ = Nozzle(cfg).emit(2000, rng)

    dirs = vel / np.linalg.norm(vel, axis=1, keepdims=True)
    axis = np.array([0.0, 0.0, -1.0])
    angles = np.arccos(np.clip(dirs @ axis, -1, 1))
    assert angles.max() <= half + 1e-6


def test_stats_are_consistent():
    cfg = SimConfig(n_droplets=1000, seed=7)
    result = Simulator(cfg).run()
    stats = analysis.summarize(result, cfg)

    assert 0.0 <= stats.landed_fraction <= 1.0
    assert stats.coverage_radius_p90 >= stats.coverage_radius_p50 >= 0.0
    assert stats.mean_radius_mm > 0.0
