"""Tests for analysis/fit_ptp_efficiency.py (coating model + fit machinery).

The simulation itself is not run here; the tests cover the model functions
and the fit on synthetic data so they stay fast.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))

import fit_ptp_efficiency as fpe  # noqa: E402


def test_uniform_shape_reduces_to_single_exponential():
    rho = np.array([0.0, 150.0, 600.0, 1100.0])
    got = fpe.model_efficiency(rho, 18.0, 400.0, np.array([1.0]))
    want = 18.0 * (1.0 - np.exp(-rho / 400.0))
    assert np.allclose(got, want)


def test_nonuniform_shape_is_below_uniform_and_saturates():
    rng = np.random.default_rng(0)
    shape = rng.lognormal(0.0, 0.8, size=5000)
    shape /= shape.mean()
    rho = np.array([100.0, 500.0, 2000.0])
    nonuni = fpe.model_efficiency(rho, 18.0, 400.0, shape)
    uni = fpe.model_efficiency(rho, 18.0, 400.0, np.array([1.0]))
    # 1 - exp(-x) is concave, so spreading the same mean thickness unevenly
    # always loses efficiency (Jensen), and both saturate at eps_sat.
    assert np.all(nonuni < uni)
    assert fpe.model_efficiency(np.array([1e6]), 18.0, 400.0, shape)[0] == pytest.approx(18.0, rel=1e-6)


def test_poisson_spot_overlap_renormalises_rho0():
    """Dried droplet spots of areal density rho_s overlapping as a Poisson
    process give exactly the single-exponential form with an effective rho0."""
    rng = np.random.default_rng(1)
    rho_s, rho0, n_mean = 60.0, 100.0, 3.0
    k = rng.poisson(n_mean, size=400_000)
    mc = np.mean(1.0 - np.exp(-k * rho_s / rho0))
    rho0_eff = fpe.poisson_effective_rho0(rho_s, rho0)
    analytic = 1.0 - np.exp(-n_mean * rho_s / rho0_eff)
    assert mc == pytest.approx(analytic, abs=2e-3)
    assert rho0_eff >= max(rho_s, rho0)


def test_fit_recovers_synthetic_parameters():
    rng = np.random.default_rng(2)
    shape = rng.lognormal(0.0, 0.5, size=2000)
    shape /= shape.mean()
    rho = np.array([150.0, 300.0, 600.0, 800.0, 1000.0, 1100.0])
    truth = (18.0, 500.0)
    eff = fpe.model_efficiency(rho, *truth, shape) + rng.normal(0.0, 0.1, size=rho.size)
    sem = np.full(rho.size, 0.1)
    fit = fpe.fit_efficiency(rho, eff, sem, shape, rho_frac_err=0.0)
    assert fit["eps_sat"] == pytest.approx(truth[0], rel=0.1)
    assert fit["rho0"] == pytest.approx(truth[1], rel=0.2)
    assert fit["ndf"] == rho.size - 2
    assert fit["chi2"] / fit["ndf"] < 3.0


def test_density_uncertainty_inflates_errors():
    rng = np.random.default_rng(3)
    shape = np.array([1.0])
    rho = np.array([150.0, 300.0, 600.0, 800.0, 1000.0, 1100.0])
    eff = fpe.model_efficiency(rho, 18.0, 500.0, shape) + rng.normal(0.0, 0.05, size=rho.size)
    sem = np.full(rho.size, 0.05)
    tight = fpe.fit_efficiency(rho, eff, sem, shape, rho_frac_err=0.0)
    loose = fpe.fit_efficiency(rho, eff, sem, shape, rho_frac_err=0.2)
    assert loose["eps_sat_err"] > tight["eps_sat_err"]
    assert loose["chi2"] < tight["chi2"]


def test_load_conf_parses_shell_assignments(tmp_path):
    conf = tmp_path / "x.conf"
    conf.write_text("# comment\nMATERIAL=toluene  # trailing\nDENSITY=\nCONE=15.0\n\nDROPLETS=1000\n")
    got = fpe.load_conf(conf)
    assert got["MATERIAL"] == "toluene"
    assert got["DENSITY"] == ""
    assert got["CONE"] == "15.0"
    assert got["DROPLETS"] == "1000"


def test_areal_density_to_thickness():
    # 123 ug/cm^2 of a 1.23 g/cm^3 solid is exactly 1 um.
    assert fpe.areal_density_to_thickness(123.0) == pytest.approx(1.0e-6)
