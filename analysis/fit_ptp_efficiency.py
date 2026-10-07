#!/usr/bin/env python3
"""Fit the IFIC pTP wavelength-shifting efficiency vs areal density with a
SpraySim coating model.

Data
----
IFIC VUV scan of p-terphenyl (pTP) sprayed from toluene onto 3 mm acrylic:
re-emitted / incident photons, QE corrected, half solid angle, 115-200 nm
(``transmittance_ptpSamples_IFIC_ok.pdf``). At a fixed wavelength (130 nm by
default) the efficiency rises with the *estimated* areal density of pTP.

Model
-----
1. **SpraySim** predicts the macroscopic dry-thickness map ``t(x, y)`` of the
   coating over the sample for the real raster toolpath and a toluene + pTP
   solution (``config/ptp_toluene.conf``). Only the *shape* of the map is used:
   it is normalised to mean 1 over the sample and scaled by the estimated
   areal density of each sample, so the non-uniformity (spray-cone profile,
   raster banding, edge fall-off) comes from the simulation and the amount of
   pTP from the data.
2. Locally, the efficiency saturates with thickness::

       eps(t) = eps_sat * (1 - exp(-t / t0)),   t0 = rho0 / rho_pTP

   ``rho0`` (ug/cm^2) lumps the VUV absorption depth of pTP with the
   droplet-scale patchiness of the dried film: a Poisson overlap of dried
   droplet spots of areal density ``rho_s`` gives *exactly* the same form with
   ``rho0 -> rho_s / (1 - exp(-rho_s / rho0))`` (see
   :func:`poisson_effective_rho0`), so the micro-scale needs no extra parameter.
3. The measured efficiency is the area average ``<eps(t)>`` over the
   illuminated region, taken as the whole sample.

Free parameters: ``eps_sat`` (saturated efficiency, %) and ``rho0`` (ug/cm^2).
Both the simulated-coating model (one per spray-cone half-angle scanned) and
the uniform-film limit (same formula, flat map) are fitted so the effect of the
simulated non-uniformity is explicit. Quoted areal densities are estimates
(dissolved pTP mass / sprayed area), so a relative density uncertainty
``--rho-frac-err`` is folded into the fit's error budget.

Usage
-----
    python analysis/fit_ptp_efficiency.py                       # defaults
    python analysis/fit_ptp_efficiency.py --cone-deg 5 10 15 25 # scan the cone
    python analysis/fit_ptp_efficiency.py --samples all         # all 14 samples
    python analysis/fit_ptp_efficiency.py --run-name ptp_130nm  # outputs in output/ptp_130nm/

Outputs go to output/<run-name>/ (output/tmp/ when no name is given).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from spraysim import (  # noqa: E402
    SimConfig, NozzleConfig, PhysicsConfig, MaterialConfig, PathConfig, Simulator, storage,
)
from spraysim.materials import (  # noqa: E402
    MATERIALS, DEFAULT_VISCOSITY, material_density, material_viscosity,
)

PTP_DENSITY = 1230.0  # kg/m^3, p-terphenyl crystal density (1.23 g/cm^3)

DEFAULT_CSV = ("/Users/razakamiandra/WORKSPACE/FD3DATA/IFIC_measurements_allSamples/"
               "transmittance_all_batches_QEcorrected.csv")
DEFAULT_CONF = ROOT / "config" / "ptp_toluene.conf"
# The six samples shown in transmittance_ptpSamples_IFIC_ok.pdf.
PLOTTED_SAMPLES = ["P067", "P070", "P048", "P054", "P060", "P064"]

# Validated categorical palette (dataviz reference instance), fixed order.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#9a9892"


# --------------------------------------------------------------------------- model

def areal_density_to_thickness(rho_ug_cm2: float | np.ndarray) -> float | np.ndarray:
    """ug/cm^2 of solid pTP -> dry thickness in metres."""
    return np.asarray(rho_ug_cm2, float) * 1e-5 / PTP_DENSITY


def model_efficiency(rho, eps_sat: float, rho0: float, shape: np.ndarray) -> np.ndarray:
    """Area-averaged efficiency (%) at areal density ``rho`` (ug/cm^2).

    ``shape`` is the normalised (mean 1) thickness map of the coating, any
    array shape; ``[1.0]`` is the uniform-film limit.
    """
    rho = np.atleast_1d(np.asarray(rho, float))
    s = np.asarray(shape, float).ravel()
    x = rho[:, None] * s[None, :] / rho0
    return eps_sat * np.mean(1.0 - np.exp(-x), axis=1)


def poisson_effective_rho0(rho_spot: float, rho0: float) -> float:
    """Effective rho0 when the film is a Poisson overlap of dried droplet spots
    of areal density ``rho_spot`` whose local response is ``1 - exp(-t/t0)``:
    <1 - exp(-k rho_s/rho0)>_{k~Poisson(n)} = 1 - exp(-n (1 - exp(-rho_s/rho0)))."""
    return rho_spot / (1.0 - math.exp(-rho_spot / rho0))


def fit_efficiency(rho, eff, sem, shape, *, rho_frac_err: float = 0.2,
                   p0=(20.0, 500.0)) -> dict:
    """Weighted least-squares fit of (eps_sat, rho0).

    The error on each point is the quoted SEM combined in quadrature with the
    model's response to a relative density error ``rho_frac_err`` (effective-
    variance method, iterated once).
    """
    from scipy.optimize import curve_fit

    rho = np.asarray(rho, float)
    eff = np.asarray(eff, float)
    sem = np.asarray(sem, float)
    s = np.asarray(shape, float).ravel()

    def f(x, eps_sat, rho0):
        return model_efficiency(x, eps_sat, rho0, s)

    sigma = sem.copy()
    popt = np.asarray(p0, float)
    for _ in range(3):
        popt, pcov = curve_fit(f, rho, eff, p0=popt, sigma=sigma, absolute_sigma=True,
                               bounds=((0.0, 1.0), (100.0, 1e6)), maxfev=20000)
        if rho_frac_err <= 0.0:
            break
        d = 1e-3 * rho
        slope = (f(rho + d, *popt) - f(rho - d, *popt)) / (2.0 * d)
        sigma = np.sqrt(sem**2 + (slope * rho_frac_err * rho) ** 2)
    resid = (eff - f(rho, *popt)) / sigma
    chi2 = float(np.sum(resid**2))
    perr = np.sqrt(np.diag(pcov))
    return {
        "eps_sat": float(popt[0]), "eps_sat_err": float(perr[0]),
        "rho0": float(popt[1]), "rho0_err": float(perr[1]),
        "t0_um": float(areal_density_to_thickness(popt[1]) * 1e6),
        "chi2": chi2, "ndf": int(rho.size - 2),
        "sigma_eff": sigma.tolist(), "residuals_pull": resid.tolist(),
    }


def fit_linear(rho, eff, sigma) -> dict:
    """Weighted straight line eps = a + b*rho — a diagnostic for the intercept.

    ``sigma`` should be the same effective errors the model fit used so the
    chi2 values are comparable.
    """
    rho = np.asarray(rho, float)
    eff = np.asarray(eff, float)
    w = 1.0 / np.asarray(sigma, float) ** 2
    A = np.vstack([np.ones_like(rho), rho]).T * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(A, eff * np.sqrt(w), rcond=None)
    chi2 = float(np.sum(w * (eff - coef[0] - coef[1] * rho) ** 2))
    return {"intercept": float(coef[0]), "slope": float(coef[1]),
            "rho_at_zero": float(-coef[0] / coef[1]) if coef[1] else float("nan"),
            "chi2": chi2, "ndf": int(rho.size - 2)}


def fit_density_offset(rho, eff, sem, *, rho_frac_err: float = 0.2) -> dict:
    """Diagnostic: uniform-film model with the estimated densities shifted by a
    common offset, eps = eps_sat (1 - exp(-(rho + d_rho)/rho0)). A d_rho far
    from 0 says the data do not extrapolate to zero efficiency at zero
    estimated density, which no coating-coverage model can reproduce."""
    from scipy.optimize import curve_fit

    rho = np.asarray(rho, float)
    eff = np.asarray(eff, float)
    sem = np.asarray(sem, float)

    def f(x, eps_sat, rho0, d_rho):
        return eps_sat * (1.0 - np.exp(-(x + d_rho) / rho0))

    sigma = sem.copy()
    popt = np.array([20.0, 500.0, 0.0])
    for _ in range(3):
        popt, pcov = curve_fit(f, rho, eff, p0=popt, sigma=sigma, absolute_sigma=True,
                               bounds=((0.0, 1.0, -5000.0), (100.0, 1e6, 5000.0)), maxfev=20000)
        if rho_frac_err <= 0.0:
            break
        d = 1e-3 * rho
        slope = (f(rho + d, *popt) - f(rho - d, *popt)) / (2.0 * d)
        sigma = np.sqrt(sem**2 + (slope * rho_frac_err * rho) ** 2)
    perr = np.sqrt(np.diag(pcov))
    chi2 = float(np.sum(((eff - f(rho, *popt)) / sigma) ** 2))
    return {"eps_sat": float(popt[0]), "eps_sat_err": float(perr[0]),
            "rho0": float(popt[1]), "rho0_err": float(perr[1]),
            "d_rho": float(popt[2]), "d_rho_err": float(perr[2]),
            "chi2": chi2, "ndf": int(rho.size - 3)}


# ---------------------------------------------------------------------- data/config

def load_data(csv_path, wavelength: float, samples: list[str] | None):
    rows = []
    with open(csv_path, newline="") as fh:
        for r in csv.DictReader(fh):
            if float(r["Wavelength (nm)"]) != wavelength:
                continue
            if samples is not None and r["Sample"] not in samples:
                continue
            rows.append({
                "sample": r["Sample"], "batch": int(r["Batch"]),
                "rho": float(r["Density (ug/cm2)"]),
                "eff": float(r["Transmission (%)"]), "sem": float(r["Transmission SEM"]),
            })
    rows.sort(key=lambda d: (d["rho"], d["sample"]))
    if not rows:
        raise SystemExit(f"no rows at {wavelength} nm in {csv_path}")
    return rows


def load_conf(path) -> dict[str, str]:
    """Parse a config/*.conf shell-style KEY=value file."""
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, val = line.split("=", 1)
        out[key.strip()] = val.strip().strip('"').strip("'")
    return out


def build_sim_config(conf: dict[str, str], *, cone_deg: float, n_droplets: int,
                     seed: int) -> SimConfig:
    name = conf.get("MATERIAL", "water")
    density = float(conf["DENSITY"]) if conf.get("DENSITY") else material_density(name)
    if conf.get("VISCOSITY"):
        viscosity = float(conf["VISCOSITY"])
    elif name in MATERIALS:
        viscosity = material_viscosity(name)
    else:
        viscosity = DEFAULT_VISCOSITY
    material = MaterialConfig(name=name, density=density, viscosity=viscosity,
                              solids_fraction=float(conf.get("SOLIDS_FRACTION", 1.0)))
    nozzle = NozzleConfig(
        half_angle=math.radians(cone_deg),
        pressure=float(conf.get("PRESSURE_BAR", 3.0)) * 1e5,
        orifice_diameter=float(conf.get("ORIFICE_MM", 0.8)) * 1e-3,
        shape=conf.get("NOZZLE_SHAPE", "full_cone"),
        speed_spread=float(conf.get("SPEED_SPREAD", 0.15)),
        distribution=conf.get("DISTRIBUTION", "lognormal"),
        mean_radius=float(conf.get("MEAN_RADIUS_MM", 0.4)) * 1e-3,
        radius_std=float(conf.get("RADIUS_STD_MM", 0.12)) * 1e-3,
    )
    gcode = conf.get("GCODE", "")
    if not gcode:
        raise SystemExit("the config must set GCODE (the raster toolpath)")
    gcode_path = Path(gcode)
    if not gcode_path.is_absolute():
        gcode_path = ROOT / gcode_path
    feed = conf.get("FEED", "")
    path = PathConfig(
        gcode=str(gcode_path),
        feed_override=float(feed) / 1000.0 / 60.0 if feed else None,
        standoff=float(conf.get("STANDOFF_MM", 150.0)) * 1e-3,
    )
    return SimConfig(
        n_droplets=n_droplets, dt=float(conf.get("DT", 1e-3)), seed=seed,
        n_trajectories=0, nozzle=nozzle, material=material, path=path,
        physics=PhysicsConfig(drag_model=conf.get("DRAG_MODEL", "clift_gauvin")),
    )


# ------------------------------------------------------------------- simulation

def simulate_shape(cfg: SimConfig, *, sample_mm: float, cell_mm: float,
                   smooth_cells: float) -> dict:
    """Run SpraySim and return the normalised thickness map over the sample.

    The sample (``sample_mm`` square) is centred on the toolpath's bounding box.
    Landed droplets are counted per cell (counts x mean droplet volume, i.e.
    the expectation over the ~1e8 real droplets rather than the sampled volume
    noise), lightly Gaussian-smoothed, then normalised to mean 1.
    """
    from scipy.ndimage import gaussian_filter

    t0 = time.time()
    result = Simulator(cfg).run()
    elapsed = time.time() - t0

    seg = result.path_segments
    xs = np.concatenate([seg[:, 0], seg[:, 2]])
    ys = np.concatenate([seg[:, 1], seg[:, 3]])
    cx, cy = 0.5 * (xs.min() + xs.max()), 0.5 * (ys.min() + ys.max())
    half = 0.5 * sample_mm * 1e-3
    cell = cell_mm * 1e-3
    n_cells = max(1, int(round(2 * half / cell)))
    x_edges = cx - half + cell * np.arange(n_cells + 1)
    y_edges = cy - half + cell * np.arange(n_cells + 1)

    pos = result.landing_positions[result.landed]
    counts, _, _ = np.histogram2d(pos[:, 0], pos[:, 1], bins=[x_edges, y_edges])
    counts = counts.T  # [iy, ix]
    in_roi = int(counts.sum())
    smooth = gaussian_filter(counts, smooth_cells, mode="nearest") if smooth_cells > 0 else counts
    mean = smooth.mean()
    shape = smooth / mean if mean > 0 else np.ones_like(smooth)

    flat = shape.ravel()
    cv = float(flat.std())
    cu = 1.0 - float(np.abs(flat - 1.0).sum()) / flat.size
    per_cell = in_roi / counts.size
    return {
        "shape": shape, "extent_mm": [1e3 * x_edges[0], 1e3 * x_edges[-1],
                                       1e3 * y_edges[0], 1e3 * y_edges[-1]],
        "cv": cv, "christiansen_cu": cu,
        "p10": float(np.percentile(flat, 10)), "p90": float(np.percentile(flat, 90)),
        "min": float(flat.min()), "max": float(flat.max()),
        "n_droplets": int(result.n), "landed_fraction": float(result.landed.mean()),
        "droplets_in_sample": in_roi, "droplets_per_cell": per_cell,
        "shot_noise_cv": (1.0 / math.sqrt(per_cell) if per_cell > 0 else float("nan")),
        "exit_speed_ms": float(result.exit_speed),
        "mean_flight_time_s": float(result.flight_times[result.landed].mean()),
        "sim_seconds": elapsed,
    }


# ------------------------------------------------------------------------- plot

def make_figure(data, others, fits, sims, best_key, wavelength, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax, axm, axh) = plt.subplots(1, 3, figsize=(15, 4.8),
                                       gridspec_kw={"width_ratios": [1.6, 1.0, 1.0]})
    fig.patch.set_facecolor("#fcfcfb")

    # --- (a) efficiency vs areal density -------------------------------------
    rho = np.array([d["rho"] for d in data])
    eff = np.array([d["eff"] for d in data])
    sem = np.array([d["sem"] for d in data])
    sig = np.array(fits[best_key]["sigma_eff"])
    grid = np.linspace(0.0, 1.15 * max(rho.max(), max((o["rho"] for o in others), default=0)), 300)

    ax.errorbar(rho, eff, yerr=sig, fmt="none", ecolor=MUTED, elinewidth=1, capsize=0, zorder=2)
    ax.errorbar(rho, eff, yerr=sem, fmt="o", color=INK, ms=6, ecolor=INK, elinewidth=1.2,
                capsize=2, zorder=4, label="fitted samples (±SEM; grey: incl. density error)")
    for d in data:
        ax.annotate(d["sample"], (d["rho"], d["eff"]), textcoords="offset points",
                    xytext=(6, -11), fontsize=8, color=INK2)
    if others:
        ax.errorbar([o["rho"] for o in others], [o["eff"] for o in others],
                    yerr=[o["sem"] for o in others], fmt="o", mfc="none", mec=MUTED,
                    ecolor=MUTED, ms=6, capsize=2, zorder=3, label="other samples (not fitted)")

    keys = [k for k in fits if k != "uniform"]
    for i, key in enumerate(keys):
        f = fits[key]
        lw = 2.4 if key == best_key else 1.6
        ax.plot(grid, model_efficiency(grid, f["eps_sat"], f["rho0"], sims[key]["shape"]),
                color=SERIES[i % len(SERIES)], lw=lw, zorder=3,
                label=(f"SpraySim coating, cone {key}°: ε$_{{sat}}$={f['eps_sat']:.1f}%, "
                       f"ρ$_0$={f['rho0']:.0f} µg/cm², χ²/ndf={f['chi2']:.1f}/{f['ndf']}"))
    fu = fits["uniform"]
    ax.plot(grid, model_efficiency(grid, fu["eps_sat"], fu["rho0"], np.array([1.0])),
            color=INK2, lw=1.6, ls="--", zorder=3,
            label=(f"uniform film: ε$_{{sat}}$={fu['eps_sat']:.1f}%, "
                   f"ρ$_0$={fu['rho0']:.0f} µg/cm², χ²/ndf={fu['chi2']:.1f}/{fu['ndf']}"))
    ax.set_xlabel("estimated pTP areal density (µg/cm²)")
    ax.set_ylabel(f"re-emitted / incident photons at {wavelength:g} nm (%)")
    ax.set_title("pTP on 3 mm acrylic — IFIC, QE corrected, half solid angle", loc="left")
    ax.set_xlim(0, grid[-1])
    ax.set_ylim(0, None)
    ax.grid(True, color="#e9e8e4", lw=0.8)
    ax.legend(fontsize=7.5, loc="lower right", frameon=False)

    # --- (b) simulated coating shape (best cone) -------------------------------
    sim = sims[best_key]
    cmap = LinearSegmentedColormap.from_list("blues", ["#fcfcfb", "#cde2fb", "#6da7ec", "#2a78d6", "#0d366b"])
    im = axm.imshow(sim["shape"], origin="lower", extent=sim["extent_mm"], cmap=cmap,
                    vmin=0, vmax=max(1.0, sim["max"]), interpolation="nearest")
    cb = fig.colorbar(im, ax=axm, fraction=0.046, pad=0.03)
    cb.set_label("local thickness / sample mean")
    cb.outline.set_edgecolor(MUTED)
    axm.set_xlabel("x on sample (mm)")
    axm.set_ylabel("y on sample (mm)")
    axm.set_title(f"SpraySim dry pTP map, cone {best_key}° (CV {sim['cv']:.2f})", loc="left")

    # --- (c) thickness distribution per cone ------------------------------------
    bins = np.linspace(0, max(s["max"] for k, s in sims.items()) * 1.02, 40)
    for i, key in enumerate(keys):
        axh.hist(sims[key]["shape"].ravel(), bins=bins, histtype="step", lw=1.8,
                 color=SERIES[i % len(SERIES)], label=f"cone {key}°: CV {sims[key]['cv']:.2f}")
    axh.axvline(1.0, color=INK2, ls="--", lw=1.2, label="uniform film")
    axh.set_xlabel("local thickness / sample mean")
    axh.set_ylabel("sample cells")
    axh.set_title("Simulated thickness spread", loc="left")
    axh.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)


# -------------------------------------------------------------------------- main

def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", default=DEFAULT_CSV, help="IFIC QE-corrected transmittance CSV")
    p.add_argument("--conf", default=str(DEFAULT_CONF), help="spray config (KEY=value)")
    p.add_argument("--wavelength", type=float, default=130.0, help="nm")
    p.add_argument("--samples", nargs="+", default=PLOTTED_SAMPLES,
                   help="samples to fit (default: the six plotted ones; 'all' = every sample)")
    p.add_argument("--cone-deg", type=float, nargs="+", default=None,
                   help="spray-cone half-angles to simulate (default: CONE from the config)")
    p.add_argument("--droplets", type=int, default=None, help="override DROPLETS from the config")
    p.add_argument("--seed", type=int, default=None, help="override SEED from the config")
    p.add_argument("--cell-mm", type=float, default=3.0, help="thickness-map cell size")
    p.add_argument("--smooth-cells", type=float, default=1.0, help="Gaussian smoothing (cells)")
    p.add_argument("--rho-frac-err", type=float, default=0.2,
                   help="relative uncertainty assumed on the estimated densities")
    p.add_argument("--run-name", default=None,
                   help="output folder name under output/ (default: 'tmp')")
    p.add_argument("--out", default="ptp_fit_130nm.png",
                   help="figure file: a bare name goes into the run folder, a path "
                        "with a directory is used as given")
    p.add_argument("--json", default=None, help="results JSON (default: alongside --out)")
    args = p.parse_args(argv)

    conf = load_conf(args.conf)
    samples = None if args.samples == ["all"] else args.samples
    data = load_data(args.csv, args.wavelength, samples)
    all_rows = load_data(args.csv, args.wavelength, None)
    fitted = {d["sample"] for d in data}
    others = [r for r in all_rows if r["sample"] not in fitted]

    cones = args.cone_deg or [float(conf.get("CONE", 25.0))]
    n_droplets = args.droplets or int(conf.get("DROPLETS") or 200000)
    seed = args.seed if args.seed is not None else int(conf.get("SEED", 42))
    sample_mm = float(conf.get("SAMPLE_MM", 74.0))

    rho = np.array([d["rho"] for d in data])
    eff = np.array([d["eff"] for d in data])
    sem = np.array([d["sem"] for d in data])

    print(f"{len(data)} samples at {args.wavelength:g} nm: "
          + ", ".join(f"{d['sample']}({d['rho']:.0f})={d['eff']:.2f}" for d in data))

    fits: dict = {"uniform": fit_efficiency(rho, eff, sem, np.array([1.0]),
                                            rho_frac_err=args.rho_frac_err)}
    sims: dict = {}
    for cone in cones:
        key = f"{cone:g}"
        cfg = build_sim_config(conf, cone_deg=cone, n_droplets=n_droplets, seed=seed)
        sim = simulate_shape(cfg, sample_mm=sample_mm, cell_mm=args.cell_mm,
                             smooth_cells=args.smooth_cells)
        sims[key] = sim
        fits[key] = fit_efficiency(rho, eff, sem, sim["shape"], rho_frac_err=args.rho_frac_err)
        print(f"cone {key}°: {sim['n_droplets']} droplets, landed {sim['landed_fraction']:.2f}, "
              f"{sim['droplets_per_cell']:.0f}/cell (shot-noise CV {sim['shot_noise_cv']:.2f}), "
              f"map CV {sim['cv']:.2f}, CU {sim['christiansen_cu']:.2f}, {sim['sim_seconds']:.0f} s")

    best_key = min((k for k in fits if k != "uniform"), key=lambda k: fits[k]["chi2"])
    linear = fit_linear(rho, eff, fits[best_key]["sigma_eff"])
    offset = fit_density_offset(rho, eff, sem, rho_frac_err=args.rho_frac_err)

    print(f"\n{'model':<18}{'eps_sat (%)':>16}{'rho0 (ug/cm2)':>20}{'t0 (um)':>10}{'chi2/ndf':>12}")
    for key, f in fits.items():
        label = "uniform film" if key == "uniform" else f"SpraySim cone {key}°"
        print(f"{label:<18}{f['eps_sat']:>9.2f} ± {f['eps_sat_err']:<5.2f}"
              f"{f['rho0']:>11.0f} ± {f['rho0_err']:<7.0f}{f['t0_um']:>8.2f}"
              f"{f['chi2']:>8.1f}/{f['ndf']}")
    print(f"diagnostics (same error model):")
    print(f"  straight line: {linear['intercept']:.2f} + {linear['slope']:.4f}*rho, "
          f"zero crossing at {linear['rho_at_zero']:.0f} ug/cm2, chi2/ndf={linear['chi2']:.1f}/{linear['ndf']}")
    print(f"  uniform film + density offset: eps_sat={offset['eps_sat']:.1f} ± {offset['eps_sat_err']:.1f} %, "
          f"rho0={offset['rho0']:.0f} ± {offset['rho0_err']:.0f}, d_rho={offset['d_rho']:.0f} ± {offset['d_rho_err']:.0f} ug/cm2, "
          f"chi2/ndf={offset['chi2']:.1f}/{offset['ndf']}")

    run_dir = storage.run_directory(args.run_name, ROOT / storage.DEFAULT_OUTPUT_ROOT)
    out_png = storage.resolve_output_path(args.out, run_dir)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    make_figure(data, others, fits, sims, best_key, args.wavelength, out_png)
    out_json = Path(args.json) if args.json else out_png.with_suffix(".json")
    payload = {
        "wavelength_nm": args.wavelength, "csv": args.csv, "conf": args.conf,
        "fitted_samples": data, "other_samples": others,
        "rho_frac_err": args.rho_frac_err, "ptp_density_kg_m3": PTP_DENSITY,
        "fits": fits, "linear": linear, "density_offset": offset, "best_cone_deg": best_key,
        "simulations": {k: {kk: vv for kk, vv in v.items() if kk != "shape"} for k, v in sims.items()},
    }
    out_json.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {out_png}\nwrote {out_json}")


if __name__ == "__main__":
    main()
