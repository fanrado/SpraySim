"""Static matplotlib plots summarising a spray run (no animation)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: render straight to file
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .config import SimConfig  # noqa: E402
from .simulator import SimResult  # noqa: E402
from . import analysis  # noqa: E402


def make_figure(result: SimResult, config: SimConfig):
    """Build a 2x2 summary figure and return the Matplotlib Figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("SpraySim — droplet spray summary", fontsize=15, fontweight="bold")

    _plot_trajectories(axes[0, 0], result, config)
    _plot_landing_scatter(axes[0, 1], result, config)
    _plot_radial_profile(axes[1, 0], result, config)
    _plot_size_histogram(axes[1, 1], result)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def _plot_trajectories(ax, result: SimResult, config: SimConfig):
    # Colour each sampled path by its droplet radius.
    sample_radii = result.radii[
        np.linspace(0, result.n - 1, len(result.trajectories), dtype=int)
    ]
    norm = plt.Normalize(sample_radii.min() * 1000, sample_radii.max() * 1000)
    cmap = plt.get_cmap("viridis")

    for traj, r in zip(result.trajectories, sample_radii):
        ax.plot(traj[:, 0], traj[:, 2], lw=0.8, alpha=0.7, color=cmap(norm(r * 1000)))

    ax.axhline(config.physics.ground_z, color="0.3", lw=1.0)
    nx, _, nz = config.nozzle.position
    ax.plot([nx], [nz], marker="v", color="crimson", ms=9, label="nozzle")
    ax.set_title("Trajectories (side view)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z / height (m)")
    ax.legend(loc="upper right", fontsize=8)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    ax.figure.colorbar(sm, ax=ax, label="radius (mm)", fraction=0.046, pad=0.04)


def _plot_landing_scatter(ax, result: SimResult, config: SimConfig):
    land = result.landing_positions[result.landed]
    speeds = result.impact_speeds[result.landed]
    sc = ax.scatter(land[:, 0], land[:, 1], c=speeds, s=6, cmap="plasma", alpha=0.6)
    nx, ny, _ = config.nozzle.position
    ax.plot([nx], [ny], marker="+", color="black", ms=12, label="spray axis")
    ax.set_title("Landing pattern (top view)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="upper right", fontsize=8)
    ax.figure.colorbar(sc, ax=ax, label="impact speed (m/s)", fraction=0.046, pad=0.04)


def _plot_radial_profile(ax, result: SimResult, config: SimConfig):
    radial = analysis.radial_distances(result, config)[result.landed]
    if radial.size == 0:
        radial = analysis.radial_distances(result, config)
    ax.hist(radial, bins=40, color="steelblue", alpha=0.85, edgecolor="white", lw=0.3)
    p50 = np.percentile(radial, 50)
    p90 = np.percentile(radial, 90)
    ax.axvline(p50, color="darkorange", ls="--", lw=1.5, label=f"p50 = {p50:.2f} m")
    ax.axvline(p90, color="crimson", ls="--", lw=1.5, label=f"p90 = {p90:.2f} m")
    ax.set_title("Radial deposition profile")
    ax.set_xlabel("distance from spray axis (m)")
    ax.set_ylabel("droplet count")
    ax.legend(fontsize=8)


def _plot_size_histogram(ax, result: SimResult):
    radii_mm = result.radii * 1000.0
    ax.hist(radii_mm, bins=40, color="seagreen", alpha=0.85, edgecolor="white", lw=0.3)
    mean_mm = radii_mm.mean()
    ax.axvline(mean_mm, color="crimson", ls="--", lw=1.5, label=f"mean = {mean_mm:.3f} mm")
    ax.set_title("Droplet size distribution")
    ax.set_xlabel("radius (mm)")
    ax.set_ylabel("droplet count")
    ax.legend(fontsize=8)


def save_figure(result: SimResult, config: SimConfig, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = make_figure(result, config)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
