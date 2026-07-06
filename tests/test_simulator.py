"""Sanity tests for the spray simulator's physics and bookkeeping."""

import math
import warnings

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
    drag,
    storage,
)
from spraysim.nozzle import Nozzle
from spraysim.simulator import SimResult

WATER = 1000.0


def _terminal_speed(radius, drag_model, height=30.0, dt=1e-3, max_time=30.0):
    """Impact speed of a droplet released from rest — i.e. its terminal speed."""
    noz = NozzleConfig(position=(0.0, 0.0, height), direction=(0, 0, -1),
                       half_angle=0.0, pressure=0.0, speed_spread=0.0,
                       distribution="normal", mean_radius=radius, radius_std=0.0)
    cfg = SimConfig(n_droplets=1, dt=dt, max_time=max_time, seed=0, nozzle=noz,
                    physics=PhysicsConfig(drag_model=drag_model))
    return Simulator(cfg).run().impact_speeds[0]


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


def test_impact_speed_energy_consistent_at_crossing():
    """Impact speed is reported at the ground crossing (P3), so in vacuum it must
    satisfy energy conservation: |v| at z=0 is sqrt(u^2 + 2 g h), independent of
    dt. Reading velocity a step late (the old behaviour) breaks this; the residual
    must also shrink with dt."""
    g, h = 9.81, 2.0

    def run(dt):
        noz = NozzleConfig(position=(0.0, 0.0, h), direction=(1, 0, 0),
                           half_angle=0.0, pressure=5.0e5, speed_spread=0.0,
                           distribution="normal", mean_radius=8e-4, radius_std=0.0)
        cfg = SimConfig(n_droplets=1, dt=dt, max_time=20.0, seed=0, nozzle=noz,
                        physics=PhysicsConfig(air_density=0.0))  # vacuum
        r = Simulator(cfg).run()
        return r.launch_speeds[0], r.impact_speeds[0]

    u, v_fine = run(1e-3)
    v_energy = math.sqrt(u * u + 2 * g * h)
    assert v_fine == pytest.approx(v_energy, rel=5e-4)
    # Finer dt lands closer to the exact energy value (O(dt) convergence).
    _, v_coarse = run(2e-3)
    assert abs(v_fine - v_energy) < abs(v_coarse - v_energy)


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


def test_timestep_convergence_is_first_order():
    """Semi-implicit Euler is O(dt): the flight-time error halves as dt halves."""
    h = 10.0

    def flight_time(dt):
        noz = NozzleConfig(position=(0.0, 0.0, h), direction=(0, 0, -1),
                           half_angle=0.0, pressure=0.0, speed_spread=0.0,
                           distribution="normal", mean_radius=4e-4, radius_std=0.0)
        cfg = SimConfig(n_droplets=1, dt=dt, max_time=30.0, seed=0, nozzle=noz,
                        physics=PhysicsConfig(drag_model="constant"))
        return Simulator(cfg).run().flight_times[0]

    ref = flight_time(1.25e-4)
    errs = [abs(flight_time(dt) - ref) for dt in (4e-3, 2e-3, 1e-3)]
    for coarse, fine in zip(errs, errs[1:]):
        assert 1.6 < coarse / fine < 2.4  # ~2x per halving


def test_smaller_droplets_travel_less_far():
    """Drag monotonicity: for a fixed launch, smaller droplets land closer."""
    ranges = []
    for r in (1e-4, 4e-4, 1.6e-3):
        noz = NozzleConfig(position=(0.0, 0.0, 2.0), direction=(1, 0, 0),
                           half_angle=0.0, pressure=3e5, speed_spread=0.0,
                           distribution="normal", mean_radius=r, radius_std=0.0)
        res = Simulator(SimConfig(n_droplets=1, dt=1e-3, max_time=20.0, seed=0,
                                  nozzle=noz)).run()
        ranges.append(math.hypot(res.landing_positions[0, 0], res.landing_positions[0, 1]))
    assert ranges[0] < ranges[1] < ranges[2]


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


def test_material_viscosity_registry_and_default():
    """Every material has a viscosity; the config default is water's viscosity."""
    assert materials.material_viscosity("water") == materials.DEFAULT_VISCOSITY
    # A viscosity is defined for every registered material.
    for name in materials.MATERIALS:
        assert materials.material_viscosity(name) > 0.0
    # Glycerin is far more viscous than water.
    assert materials.material_viscosity("glycerin") > materials.material_viscosity("water")
    # Unspecified material viscosity defaults to water's.
    assert MaterialConfig().viscosity == materials.material_viscosity("water")
    with pytest.raises(KeyError):
        materials.material_viscosity("unobtanium")


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
    """mean_cubed_radius matches the sampled E[r^3], including a wide (heavily
    clipped) normal where the naive m^3+3ms^2 would drift ~8%."""
    rng = np.random.default_rng(0)
    cases = [("normal", 4e-4, 1e-4), ("normal", 2e-4, 3e-4), ("lognormal", 4e-4, 1e-4)]
    for dist, m, s in cases:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # the wide normal warns; not under test here
            nozzle = Nozzle(NozzleConfig(distribution=dist, mean_radius=m, radius_std=s), WATER)
            sampled = np.mean(nozzle.sample_radii(500_000, rng) ** 3)
        assert nozzle.mean_cubed_radius() == pytest.approx(sampled, rel=2e-2)


def test_wide_normal_distribution_warns():
    """A wide 'normal' spread (many radii clipped at 0) warns; narrow normal and
    lognormal stay silent."""
    with pytest.warns(UserWarning, match="clips"):
        Nozzle(NozzleConfig(distribution="normal", mean_radius=2e-4, radius_std=2e-4), WATER)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning here fails the test
        Nozzle(NozzleConfig(distribution="normal", mean_radius=4e-4, radius_std=1e-4), WATER)
        Nozzle(NozzleConfig(distribution="lognormal", mean_radius=2e-4, radius_std=3e-4), WATER)


def test_npz_round_trip_preserves_result_and_config(tmp_path):
    """Saving to .npz and loading back reproduces arrays, trajectories, config."""
    cfg = SimConfig(n_droplets=300, seed=11,
                    material=MaterialConfig("ethanol", 789.0, 1.2e-3, solids_fraction=0.3))
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
    assert loaded_cfg.material.viscosity == pytest.approx(cfg.material.viscosity)
    assert loaded_cfg.material.solids_fraction == pytest.approx(cfg.material.solids_fraction)
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


def test_drag_coefficient_limits():
    """Clift-Gauvin reduces to Stokes (24/Re) at low Re and ~Newton at high Re."""
    lo = float(drag.clift_gauvin(np.array(0.01)))
    assert lo == pytest.approx(24.0 / 0.01, rel=0.02)  # Stokes limit
    hi = float(drag.clift_gauvin(np.array(1e5)))
    assert 0.35 < hi < 0.55                              # Newton plateau
    # "constant" model ignores Re and returns the given coefficient.
    assert float(drag.drag_coefficient(np.array(123.0), "constant", 0.47)) == 0.47
    with pytest.raises(ValueError):
        drag.drag_coefficient(np.array(1.0), "no_such_model")


def test_constant_model_matches_legacy_terminal_velocity():
    """The 'constant' model reproduces the pre-P1 terminal velocity v_t=sqrt(g/k)."""
    r = 8e-4
    k = 3.0 * 1.225 * 0.47 / (8.0 * WATER * r)
    v_expected = math.sqrt(9.81 / k)
    assert _terminal_speed(r, "constant") == pytest.approx(v_expected, rel=1e-3)


def test_reynolds_drag_matches_fixed_point_and_slows_small_drops():
    """clift_gauvin terminal speed matches a Cd(Re) fixed-point solve, and is
    below the constant-Cd value (fine droplets are slowed more)."""
    r = 8e-4
    # Fixed-point solve of g = k(Cd(Re)) v^2 with the same correlation.
    v = 5.0
    for _ in range(200):
        Re = 1.225 * v * 2 * r / drag.DEFAULT_AIR_VISCOSITY
        cd = float(drag.clift_gauvin(np.array(Re)))
        v = math.sqrt(9.81 / (3.0 * 1.225 * cd / (8.0 * WATER * r)))
    assert _terminal_speed(r, "clift_gauvin") == pytest.approx(v, rel=1e-2)
    # More drag than the constant model for a small droplet.
    assert _terminal_speed(2e-4, "clift_gauvin") < _terminal_speed(2e-4, "constant")


def test_drag_fields_round_trip_and_old_archive_fallback(tmp_path):
    """air_viscosity + drag_model round-trip; archives lacking them load as
    the legacy constant-Cd model."""
    cfg = SimConfig(n_droplets=50, seed=1,
                    physics=PhysicsConfig(drag_model="clift_gauvin",
                                          air_viscosity=1.9e-5))
    result = Simulator(cfg).run()
    path = storage.save_result(result, cfg, tmp_path / "run.npz")
    _, loaded = storage.load_result(path)
    assert loaded.physics.drag_model == "clift_gauvin"
    assert loaded.physics.air_viscosity == pytest.approx(1.9e-5)

    # Simulate a pre-P1 archive by stripping the new keys.
    data = {k: v for k, v in np.load(path).items()
            if k not in ("cfg_drag_model", "cfg_air_viscosity")}
    old = tmp_path / "legacy.npz"
    np.savez(old, **data)
    _, legacy = storage.load_result(old)
    assert legacy.physics.drag_model == "constant"
    assert legacy.physics.air_viscosity == pytest.approx(drag.DEFAULT_AIR_VISCOSITY)


def test_stats_are_consistent():
    cfg = SimConfig(n_droplets=1000, seed=7)
    result = Simulator(cfg).run()
    stats = analysis.summarize(result, cfg)

    assert 0.0 <= stats.landed_fraction <= 1.0
    assert stats.coverage_radius_p90 >= stats.coverage_radius_p50 >= 0.0
    assert stats.mean_radius_mm > 0.0


def _landed_result(positions, radii):
    """A minimal SimResult with the given landings (all marked landed)."""
    n = len(radii)
    return SimResult(
        landing_positions=np.asarray(positions, dtype=float),
        flight_times=np.zeros(n),
        impact_speeds=np.zeros(n),
        radii=np.asarray(radii, dtype=float),
        launch_speeds=np.zeros(n),
        trajectories=[],
        landed=np.ones(n, dtype=bool),
    )


def test_deposition_map_conserves_volume_and_thickness():
    """A single cell over a known area gives thickness = deposited volume / area,
    and the gridded wet volume equals the summed droplet volume."""
    rng = np.random.default_rng(0)
    n, r0, L = 1000, 2e-4, 0.5
    xy = rng.uniform(0.0, L, size=(n, 2))
    pos = np.column_stack([xy, np.zeros(n)])
    result = _landed_result(pos, np.full(n, r0))
    cfg = SimConfig(material=MaterialConfig(solids_fraction=1.0))

    field = analysis.deposition_map(result, cfg, cell_size=L, extent=(0, L, 0, L))
    droplet_vol = n * (4.0 / 3.0) * np.pi * r0 ** 3
    assert field.thickness.shape == (1, 1)
    assert field.total_wet_volume() == pytest.approx(droplet_vol, rel=1e-9)
    assert field.thickness[0, 0] == pytest.approx(droplet_vol / L ** 2, rel=1e-9)


def test_solids_fraction_scales_dry_thickness_only():
    """Dry thickness scales with solids_fraction; wet volume is unchanged."""
    rng = np.random.default_rng(1)
    n, r0, L, sf = 500, 3e-4, 0.4, 0.25
    pos = np.column_stack([rng.uniform(0, L, (n, 2)), np.zeros(n)])
    result = _landed_result(pos, np.full(n, r0))
    droplet_vol = n * (4.0 / 3.0) * np.pi * r0 ** 3

    dry = analysis.deposition_map(
        result, SimConfig(material=MaterialConfig(solids_fraction=sf)),
        cell_size=L, extent=(0, L, 0, L))
    assert dry.thickness[0, 0] == pytest.approx(sf * droplet_vol / L ** 2, rel=1e-9)
    assert dry.total_dry_volume() == pytest.approx(sf * droplet_vol, rel=1e-9)
    assert dry.total_wet_volume() == pytest.approx(droplet_vol, rel=1e-9)


def test_deposition_map_total_matches_real_run():
    """On a real run, the auto-extent grid captures the full sprayed volume."""
    cfg = SimConfig(n_droplets=800, seed=4)
    result = Simulator(cfg).run()
    field = analysis.deposition_map(result, cfg)
    landed = result.landed
    expected = np.sum((4.0 / 3.0) * np.pi * result.radii[landed] ** 3)
    assert field.total_wet_volume() == pytest.approx(expected, rel=1e-6)


def _uniform_field(value, ny=10, nx=10, cell=1e-2):
    t = np.full((ny, nx), float(value))
    x_edges = np.arange(nx + 1) * cell
    y_edges = np.arange(ny + 1) * cell
    return analysis.DepositionField(x_edges, y_edges, t, cell, 1.0)


def test_uniformity_of_uniform_field_is_perfect():
    u = analysis.uniformity(_uniform_field(5e-6))
    assert u.cv == pytest.approx(0.0, abs=1e-12)
    assert u.christiansen_cu == pytest.approx(1.0, abs=1e-12)
    assert u.coverage_fraction == pytest.approx(1.0)
    assert u.mean_thickness == pytest.approx(5e-6)


def test_uniformity_single_spot_is_nonuniform():
    """A single-spot spray has strong radial fall-off: CV high, CU well below 1."""
    cfg = SimConfig(n_droplets=3000, seed=7)
    result = Simulator(cfg).run()
    u = analysis.uniformity(analysis.deposition_map(result, cfg))
    assert u.n_cells > 0
    assert 0.0 <= u.coverage_fraction <= 1.0
    assert u.cv > 0.1
    assert u.christiansen_cu < 0.9


def test_uniformity_rectangle_roi_penalizes_gaps():
    """Over a target rectangle, an uncoated half drops coverage to ~0.5; the
    default (wetted) ROI sees only the coated half and reads uniform."""
    ny, nx, cell = 10, 10, 1e-2
    t = np.zeros((ny, nx))
    t[:, :5] = 1e-6  # coat the left half only
    x_edges = np.arange(nx + 1) * cell
    y_edges = np.arange(ny + 1) * cell
    field = analysis.DepositionField(x_edges, y_edges, t, cell, 1.0)

    full = analysis.uniformity(field, roi=(0.0, nx * cell, 0.0, ny * cell))
    assert full.coverage_fraction == pytest.approx(0.5, abs=0.05)
    assert full.christiansen_cu < 1.0

    wetted = analysis.uniformity(field)  # default ROI = non-zero cells
    assert wetted.christiansen_cu == pytest.approx(1.0, abs=1e-12)
    assert wetted.coverage_fraction == pytest.approx(1.0)


def test_plot_deposition_smoke():
    """plot_deposition renders a heatmap figure without error."""
    import matplotlib.pyplot as plt
    from spraysim import plots

    fig = plots.plot_deposition(_uniform_field(3e-6))
    assert fig is not None
    plt.close(fig)
