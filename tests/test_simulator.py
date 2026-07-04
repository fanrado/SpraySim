"""Sanity tests for the spray simulator's physics and bookkeeping."""

import math

import numpy as np
import pytest

from spraysim import (
    SimConfig,
    NozzleConfig,
    PhysicsConfig,
    MaterialConfig,
    Simulator,
    analysis,
    hydraulics,
    materials,
    storage,
)
from spraysim.nozzle import Nozzle

WATER = 1000.0


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
    # pressure=0 => exit speed 0 => a droplet released from rest.
    nozzle = NozzleConfig(position=(0.0, 0.0, h), half_angle=0.0, pressure=0.0,
                          speed_spread=0.0)
    cfg = SimConfig(n_droplets=1, dt=1e-4, nozzle=nozzle, physics=phys, seed=0)

    result = Simulator(cfg).run()
    expected = math.sqrt(2 * h / phys.gravity)
    assert result.exit_speed == pytest.approx(0.0)
    assert result.flight_times[0] == pytest.approx(expected, rel=1e-2)


def test_cone_directions_respect_half_angle():
    half = math.radians(20.0)
    cfg = NozzleConfig(direction=(0, 0, -1), half_angle=half, speed_spread=0.0)
    rng = np.random.default_rng(0)
    _, vel, _ = Nozzle(cfg, WATER).emit(2000, rng)

    dirs = vel / np.linalg.norm(vel, axis=1, keepdims=True)
    axis = np.array([0.0, 0.0, -1.0])
    angles = np.arccos(np.clip(dirs @ axis, -1, 1))
    assert angles.max() <= half + 1e-6


def test_exit_speed_scales_with_sqrt_pressure():
    """Torricelli: doubling pressure raises exit speed by sqrt(2)."""
    v1 = hydraulics.exit_speed(2.0e5, "full_cone", WATER)
    v2 = hydraulics.exit_speed(4.0e5, "full_cone", WATER)
    assert v2 / v1 == pytest.approx(math.sqrt(2.0), rel=1e-9)


def test_droplet_count_derivation_matches_volume_balance():
    """Derived count ~= sprayed volume / mean droplet volume."""
    cfg = SimConfig(spray_duration=0.2, seed=3)  # n_droplets=None => derive
    result = Simulator(cfg).run()

    nozzle = Nozzle(cfg.nozzle, cfg.material.density)
    expected = nozzle.flow_rate * cfg.spray_duration / nozzle.mean_droplet_volume()
    assert result.n == pytest.approx(round(expected))
    assert not result.droplets_capped
    assert result.flow_rate > 0.0


def test_more_pressure_yields_more_droplets():
    low = NozzleConfig(pressure=2.0e5)
    high = NozzleConfig(pressure=6.0e5)
    n_low = Simulator(SimConfig(nozzle=low)).run().n
    n_high = Simulator(SimConfig(nozzle=high)).run().n
    assert n_high > n_low


def test_material_density_registry():
    """Named materials resolve to a density; unknown names raise clearly."""
    assert materials.material_density("water") == 1000.0
    assert materials.material_density("ethanol") < materials.material_density("water")
    with pytest.raises(KeyError):
        materials.material_density("unobtanium")


def test_denser_liquid_lowers_exit_speed():
    """Torricelli: v = C_v sqrt(2 dP / rho), so exit speed scales as 1/sqrt(rho)."""
    light = hydraulics.exit_speed(3.0e5, "full_cone", 800.0)
    heavy = hydraulics.exit_speed(3.0e5, "full_cone", 1200.0)
    assert heavy < light
    assert light / heavy == pytest.approx(math.sqrt(1200.0 / 800.0), rel=1e-9)


def test_material_changes_simulation_outcome():
    """Switching the sprayed liquid changes droplet dynamics for a fixed nozzle."""
    base = SimConfig(n_droplets=400, seed=2, material=MaterialConfig("water", 1000.0))
    dense = SimConfig(n_droplets=400, seed=2, material=MaterialConfig("glycerin", 1260.0))
    r_water = Simulator(base).run()
    r_glyc = Simulator(dense).run()
    # Denser liquid -> lower Torricelli exit speed for the same pressure.
    assert r_glyc.exit_speed < r_water.exit_speed
    # The landing pattern is not identical between materials.
    assert not np.allclose(r_glyc.landing_positions, r_water.landing_positions)


def test_normal_and_lognormal_moments():
    """mean_cubed_radius matches the sampled E[r^3] for both distributions."""
    rng = np.random.default_rng(0)
    for dist in ("normal", "lognormal"):
        cfg = NozzleConfig(distribution=dist, mean_radius=4e-4, radius_std=1e-4)
        nozzle = Nozzle(cfg, WATER)
        sampled = np.mean(nozzle.sample_radii(200_000, rng) ** 3)
        assert nozzle.mean_cubed_radius() == pytest.approx(sampled, rel=2e-2)


def test_npz_round_trip_preserves_result_and_config(tmp_path):
    """Saving to .npz and loading back reproduces arrays, trajectories, config."""
    cfg = SimConfig(n_droplets=300, seed=11,
                    material=MaterialConfig("ethanol", 789.0))
    result = Simulator(cfg).run()

    path = storage.save_result(result, cfg, tmp_path / "run.npz")
    assert path.exists()

    loaded, loaded_cfg = storage.load_result(path)

    # Per-droplet arrays survive exactly.
    assert np.array_equal(loaded.landing_positions, result.landing_positions)
    assert np.array_equal(loaded.radii, result.radii)
    assert np.array_equal(loaded.landed, result.landed)
    assert loaded.exit_speed == pytest.approx(result.exit_speed)
    assert loaded.droplets_capped == result.droplets_capped

    # Variable-length trajectories split back to the same shapes and values.
    assert len(loaded.trajectories) == len(result.trajectories)
    for a, b in zip(loaded.trajectories, result.trajectories):
        assert np.array_equal(a, b)

    # Config round-trips, so re-running the loaded config reproduces the result.
    assert loaded_cfg.n_droplets == cfg.n_droplets
    assert loaded_cfg.seed == cfg.seed
    assert loaded_cfg.nozzle.shape == cfg.nozzle.shape
    assert loaded_cfg.nozzle.distribution == cfg.nozzle.distribution
    assert loaded_cfg.nozzle.pressure == pytest.approx(cfg.nozzle.pressure)
    assert loaded_cfg.material.name == cfg.material.name
    assert loaded_cfg.material.density == pytest.approx(cfg.material.density)
    rerun = Simulator(loaded_cfg).run()
    assert np.array_equal(rerun.landing_positions, result.landing_positions)


def test_npz_round_trip_with_derived_count(tmp_path):
    """A derived-count run (n_droplets=None) round-trips the None sentinel."""
    cfg = SimConfig(spray_duration=0.05, seed=5)  # n_droplets=None => derived
    result = Simulator(cfg).run()

    path = storage.save_result(result, cfg, tmp_path / "derived.npz")
    _, loaded_cfg = storage.load_result(path)

    assert loaded_cfg.n_droplets is None
    assert loaded_cfg.spray_duration == pytest.approx(cfg.spray_duration)


def test_stats_are_consistent():
    cfg = SimConfig(n_droplets=1000, seed=7)
    result = Simulator(cfg).run()
    stats = analysis.summarize(result, cfg)

    assert 0.0 <= stats.landed_fraction <= 1.0
    assert stats.coverage_radius_p90 >= stats.coverage_radius_p50 >= 0.0
    assert stats.mean_radius_mm > 0.0
