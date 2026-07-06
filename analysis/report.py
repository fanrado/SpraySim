#!/usr/bin/env python3
"""Build a PDF report from saved SpraySim runs (``output/*.npz``).

Each ``.npz`` archive written by a simulation run (see ``spraysim.storage``) is
fully self-describing — it carries the per-droplet arrays *and* the config that
produced them — so this script can regenerate every plot and statistic offline,
without re-running the simulation.

Usage
-----
    python analysis/report.py                       # all output/*.npz -> output/spray_report.pdf
    python analysis/report.py output/big_drops.npz  # one run
    python analysis/report.py a.npz b.npz --out report.pdf
    python analysis/report.py --glob "output/*.npz" --out output/spray_report.pdf

The PDF contains, per run, the standard 2x2 summary figure, an extra-analysis
page (radial coverage CDF and size-vs-range scatter) and a config/stats page;
plus a cover page and, when several runs are given, a comparison table.
"""

from __future__ import annotations

import argparse
import glob as globmod
import sys
from datetime import datetime
from pathlib import Path

# Allow "python analysis/report.py" from the repo root without installing.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from spraysim import analysis, plots, storage  # noqa: E402


def load_runs(paths: list[Path]):
    """Load each .npz into (name, result, config, stats).

    Archives that cannot be read (e.g. written by an older, incompatible format)
    are skipped with a warning so one bad file does not sink the whole report.
    """
    runs = []
    for p in paths:
        try:
            result, config = storage.load_result(p)
        except (KeyError, ValueError, OSError) as exc:
            print(f"warning: skipping {p} (not a readable SpraySim archive: {exc})",
                  file=sys.stderr)
            continue
        stats = analysis.summarize(result, config)
        runs.append((p.stem, result, config, stats))
    if not runs:
        raise SystemExit("error: none of the given .npz files were readable runs.")
    return runs


# --------------------------------------------------------------------------- #
# Pages
# --------------------------------------------------------------------------- #

def _text_page(pdf: PdfPages, title: str, lines: list[str]) -> None:
    """Render a monospaced text page (config / stats dumps)."""
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    fig.text(0.5, 0.95, title, ha="center", va="top", fontsize=16, fontweight="bold")
    fig.text(0.08, 0.90, "\n".join(lines), ha="left", va="top",
             fontsize=10, family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def cover_page(pdf: PdfPages, runs) -> None:
    lines = [
        f"Generated: {datetime.now():%Y-%m-%d %H:%M}",
        f"Runs in this report: {len(runs)}",
        "",
    ]
    for name, result, config, _ in runs:
        lines.append(f"  - {name}: {config.material.name}, {result.n} droplets, "
                     f"{config.nozzle.pressure / 1e5:g} bar")
    _text_page(pdf, "SpraySim — Analysis Report", lines)


def comparison_page(pdf: PdfPages, runs) -> None:
    """A table comparing key metrics across runs (only when >1 run)."""
    cols = ["material", "droplets", "exit v\n(m/s)", "flow\n(mL/s)",
            "mean t\n(s)", "p90 r\n(m)", "mean rad\n(mm)", "CU"]
    rows, cells = [], []
    for name, result, config, stats in runs:
        cu = analysis.uniformity(analysis.deposition_map(result, config)).christiansen_cu
        rows.append(name)
        cells.append([
            config.material.name,
            f"{result.n}",
            f"{result.exit_speed:.2f}",
            f"{result.flow_rate * 1e6:.2f}",
            f"{stats.mean_flight_time:.3f}",
            f"{stats.coverage_radius_p90:.3f}",
            f"{stats.mean_radius_mm:.3f}",
            f"{cu:.2f}",
        ])

    fig, ax = plt.subplots(figsize=(11.69, 8.27))  # A4 landscape
    ax.axis("off")
    ax.set_title("Run comparison", fontsize=15, fontweight="bold", pad=20)
    table = ax.table(cellText=cells, rowLabels=rows, colLabels=cols,
                     loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.8)
    pdf.savefig(fig)
    plt.close(fig)


def _extra_analysis_page(pdf: PdfPages, name, result, config) -> None:
    """Radial coverage CDF and droplet-size-vs-range scatter."""
    radial = analysis.radial_distances(result, config)
    landed = result.landed
    r_land = radial[landed] if landed.any() else radial
    radii_mm = result.radii * 1000.0
    radii_land = radii_mm[landed] if landed.any() else radii_mm

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"{name} — coverage & size sorting", fontsize=14, fontweight="bold")

    # Cumulative coverage: fraction of droplets landing within radius r.
    order = np.sort(r_land)
    cdf = np.arange(1, order.size + 1) / order.size
    ax1.plot(order, cdf, color="steelblue", lw=1.8)
    for q, c in ((0.5, "darkorange"), (0.9, "crimson")):
        rq = float(np.percentile(r_land, q * 100))
        ax1.axvline(rq, color=c, ls="--", lw=1.2, label=f"p{int(q*100)} = {rq:.2f} m")
    ax1.set_title("Radial coverage (CDF)")
    ax1.set_xlabel("distance from spray axis (m)")
    ax1.set_ylabel("fraction of droplets within r")
    ax1.set_ylim(0, 1)
    ax1.legend(fontsize=8)

    # Size sorting: do bigger droplets fly further?
    sc = ax2.scatter(radii_land, r_land, s=6, alpha=0.5,
                     c=result.impact_speeds[landed] if landed.any() else result.impact_speeds,
                     cmap="plasma")
    ax2.set_title("Landing distance vs droplet size")
    ax2.set_xlabel("droplet radius (mm)")
    ax2.set_ylabel("distance from spray axis (m)")
    fig.colorbar(sc, ax=ax2, label="impact speed (m/s)", fraction=0.046, pad=0.04)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def _config_stats_page(pdf: PdfPages, name, result, config, stats) -> None:
    noz, mat, phys = config.nozzle, config.material, config.physics
    lines = [
        "MATERIAL",
        f"  name          {mat.name}",
        f"  density       {mat.density:g} kg/m^3",
        f"  viscosity     {mat.viscosity:g} Pa*s",
        "",
        "NOZZLE / HYDRAULICS",
        f"  pressure      {noz.pressure / 1e5:g} bar",
        f"  orifice       {noz.orifice_diameter * 1e3:g} mm",
        f"  shape         {noz.shape}",
        f"  cone half-ang {np.degrees(noz.half_angle):g} deg",
        f"  height        {noz.position[2]:g} m",
        f"  distribution  {noz.distribution} "
        f"(mean {noz.mean_radius * 1e3:g} mm, std {noz.radius_std * 1e3:g} mm)",
        f"  spray time    {config.spray_duration:g} s",
        "",
        "DERIVED",
        f"  exit speed    {result.exit_speed:.3f} m/s",
        f"  flow rate     {result.flow_rate * 1e6:.3f} mL/s",
        f"  droplets      {result.n}"
        + ("  (CAPPED)" if result.droplets_capped else ""),
        "",
        "STATISTICS",
    ]
    for k, v in stats.as_dict().items():
        lines.append(f"  {k:<24} {_fmt(v)}")
    _text_page(pdf, f"{name} — configuration & statistics", lines)


def _fmt(v):
    """Compactly format a stat value (round floats, keep ints/strings)."""
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (tuple, list)):
        return "(" + ", ".join(_fmt(x) for x in v) + ")"
    return str(v)


def _deposition_page(pdf: PdfPages, name, result, config) -> None:
    """Dry film-thickness heatmap + per-cell thickness histogram with CU/CV."""
    field = analysis.deposition_map(result, config)
    u = analysis.uniformity(field)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"{name} — deposition & uniformity", fontsize=14, fontweight="bold")
    plots.plot_deposition(field, ax=ax1)

    vals = field.thickness[field.nonzero_mask()] * 1e6
    if vals.size:
        ax2.hist(vals, bins=40, color="steelblue", alpha=0.85, edgecolor="white", lw=0.3)
        ax2.axvline(u.mean_thickness * 1e6, color="crimson", ls="--", lw=1.5,
                    label=f"mean = {u.mean_thickness * 1e6:.3f} µm")
        ax2.legend(fontsize=8)
    ax2.set_title(f"Cell thickness — CU={u.christiansen_cu:.2f}, "
                  f"CV={u.cv:.2f}, coverage={u.coverage_fraction:.0%}")
    ax2.set_xlabel("dry thickness (µm)")
    ax2.set_ylabel("cell count")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def run_pages(pdf: PdfPages, name, result, config, stats) -> None:
    fig = plots.make_figure(result, config)
    fig.suptitle(f"SpraySim — {name}", fontsize=15, fontweight="bold")
    pdf.savefig(fig)
    plt.close(fig)
    _extra_analysis_page(pdf, name, result, config)
    _deposition_page(pdf, name, result, config)
    _config_stats_page(pdf, name, result, config, stats)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def build_report(paths: list[Path], out: Path) -> Path:
    runs = load_runs(paths)
    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        cover_page(pdf, runs)
        if len(runs) > 1:
            comparison_page(pdf, runs)
        for name, result, config, stats in runs:
            run_pages(pdf, name, result, config, stats)
    return out


def resolve_paths(args: argparse.Namespace) -> list[Path]:
    if args.inputs:
        paths = [Path(p) for p in args.inputs]
    else:
        paths = sorted(Path(p) for p in globmod.glob(args.glob))
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise SystemExit(f"error: file(s) not found: {', '.join(map(str, missing))}")
    if not paths:
        raise SystemExit(
            f"error: no .npz files matched {args.glob!r}. Run a simulation first "
            "(e.g. ./main.sh) to produce output/*.npz."
        )
    return paths


def main() -> None:
    p = argparse.ArgumentParser(description="Generate a PDF report from output/*.npz runs.")
    p.add_argument("inputs", nargs="*", help="specific .npz files (default: --glob)")
    p.add_argument("--glob", default="output/*.npz",
                   help="glob for input archives when none are named (default: output/*.npz)")
    p.add_argument("--out", type=Path, default=Path("output/spray_report.pdf"),
                   help="output PDF path (default: output/spray_report.pdf)")
    args = p.parse_args()

    paths = resolve_paths(args)
    out = build_report(paths, args.out)
    print(f"Report written to {out.resolve()} ({len(paths)} run(s))")


if __name__ == "__main__":
    main()
