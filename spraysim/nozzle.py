"""Droplet emitter: samples initial positions, velocities and radii.

The exit speed is derived from the nozzle pressure and shape via
:mod:`spraysim.hydraulics`; only the droplet-size distribution is sampled here.
"""

from __future__ import annotations

import math
import warnings

import numpy as np

from .config import NozzleConfig
from . import hydraulics

_MIN_RADIUS = 1.0e-9  # m, floor to keep radii strictly positive after clipping
# Warn when the normal distribution clips more than this fraction of radii at 0.
_CLIP_WARN_FRACTION = 0.01


def _orthonormal_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u, v, w) with w == normalised axis and u, v spanning its plane."""
    w = axis / np.linalg.norm(axis)
    # Pick a helper vector that is not parallel to w.
    helper = np.array([1.0, 0.0, 0.0]) if abs(w[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, w)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)
    return u, v, w


def _allocate_counts(n: int, weights: np.ndarray) -> np.ndarray:
    """Split ``n`` items across bins by ``weights``, summing exactly to ``n``.

    Largest-remainder rounding: floor each share, then hand the leftover items to
    the bins with the biggest fractional parts.
    """
    exact = n * weights
    base = np.floor(exact).astype(int)
    remainder = int(n - base.sum())
    if remainder > 0:
        order = np.argsort(exact - base)[::-1][:remainder]
        base[order] += 1
    return base


def _lognormal_params(mean: float, std: float) -> tuple[float, float]:
    """Underlying-normal (mu, sigma) of a log-normal with the given linear mean/std."""
    sigma2 = math.log(1.0 + (std / mean) ** 2)
    sigma = math.sqrt(sigma2)
    mu = math.log(mean) - 0.5 * sigma2
    return mu, sigma


class Nozzle:
    """Emits droplets into a cone around a spray axis."""

    def __init__(self, config: NozzleConfig, liquid_density: float):
        self.config = config
        self.liquid_density = liquid_density
        self._u, self._v, self._w = _orthonormal_basis(np.asarray(config.direction, float))

        # Derived hydraulics (constant for the run).
        self.exit_speed = hydraulics.exit_speed(config.pressure, config.shape, liquid_density)
        self.flow_rate = hydraulics.flow_rate(
            config.pressure, config.orifice_diameter, config.shape, liquid_density
        )

        self._warn_if_normal_clips()

    def _warn_if_normal_clips(self) -> None:
        """Warn when a wide 'normal' distribution clips many radii at ~0.

        E[r^3] (and thus the droplet count) is corrected for the clipping, but a
        large clipped fraction still means many unphysical near-zero droplets;
        'lognormal' avoids negative radii entirely.
        """
        cfg = self.config
        if cfg.distribution != "normal" or cfg.radius_std <= 0.0:
            return
        # Fraction of the Gaussian below zero: Phi(-mean/std).
        clip_frac = 0.5 * math.erfc((cfg.mean_radius / cfg.radius_std) / math.sqrt(2.0))
        if clip_frac > _CLIP_WARN_FRACTION:
            warnings.warn(
                f"normal droplet-size distribution clips {clip_frac:.1%} of radii "
                f"at 0 (mean_radius/radius_std="
                f"{cfg.mean_radius / cfg.radius_std:.2f}); consider "
                "distribution='lognormal' for a wide spread.",
                stacklevel=3,
            )

    # -- droplet size distribution -------------------------------------------

    def sample_radii(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Draw ``n`` droplet radii (m) from the configured distribution."""
        cfg = self.config
        if cfg.distribution == "normal":
            radii = rng.normal(cfg.mean_radius, cfg.radius_std, size=n)
            return np.clip(radii, _MIN_RADIUS, None)
        if cfg.distribution == "lognormal":
            mu, sigma = _lognormal_params(cfg.mean_radius, cfg.radius_std)
            return rng.lognormal(mean=mu, sigma=sigma, size=n)
        raise ValueError(
            f"Unknown distribution {cfg.distribution!r}; use 'normal' or 'lognormal'."
        )

    def mean_cubed_radius(self) -> float:
        """Analytic E[r^3] of the size distribution (drives droplet volume)."""
        cfg = self.config
        m, s = cfg.mean_radius, cfg.radius_std
        if cfg.distribution == "normal":
            # The sampler clips radii at ~0, so the relevant moment is the
            # partial (left-truncated at 0) third moment E[max(r,0)^3], which
            # matches the clipped samples for any spread. The full m^3+3 m s^2
            # drifts high by the neglected negative tail once s is not << m.
            # For X ~ N(m, s^2):
            #   Phi(m/s) * (m^3 + 3 m s^2) + phi(m/s) * s * (m^2 + 2 s^2)
            if s <= 0.0:
                return m**3
            z = m / s
            Phi = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
            phi = math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)
            return Phi * (m**3 + 3.0 * m * s**2) + phi * s * (m**2 + 2.0 * s**2)
        if cfg.distribution == "lognormal":
            mu, sigma = _lognormal_params(m, s)
            return math.exp(3.0 * mu + 4.5 * sigma**2)
        raise ValueError(
            f"Unknown distribution {cfg.distribution!r}; use 'normal' or 'lognormal'."
        )

    def mean_droplet_volume(self) -> float:
        """Mean volume (m^3) of a droplet: (4/3) pi E[r^3]."""
        return (4.0 / 3.0) * math.pi * self.mean_cubed_radius()

    # -- emission ------------------------------------------------------------

    def _emit_core(self, positions: np.ndarray, rng: np.random.Generator,
                   carriage_velocity: np.ndarray | None = None):
        """Sample droplet velocities/radii for the given launch ``positions``.

        ``carriage_velocity`` ((3,) or (n, 3), m/s) is added to each droplet's
        launch velocity — a moving nozzle throws droplets along its travel.
        """
        cfg = self.config
        n = positions.shape[0]

        # Directions: uniform over the cone's solid angle.
        # cos(theta) uniform in [cos(half_angle), 1] gives an even areal spread.
        cos_max = np.cos(cfg.half_angle)
        cos_theta = rng.uniform(cos_max, 1.0, size=n)
        sin_theta = np.sqrt(np.clip(1.0 - cos_theta**2, 0.0, 1.0))
        phi = rng.uniform(0.0, 2.0 * np.pi, size=n)

        # Direction in nozzle frame -> world frame via the orthonormal basis.
        dirs = (
            (sin_theta * np.cos(phi))[:, None] * self._u
            + (sin_theta * np.sin(phi))[:, None] * self._v
            + cos_theta[:, None] * self._w
        )

        # Speeds: positive, spread around the pressure-derived exit speed.
        speeds = rng.normal(self.exit_speed, self.exit_speed * cfg.speed_spread, size=n)
        speeds = np.clip(speeds, 0.0, None)
        velocities = dirs * speeds[:, None]
        if carriage_velocity is not None:
            velocities = velocities + carriage_velocity

        radii = self.sample_radii(n, rng)
        return positions, velocities, radii

    def emit(self, n: int, rng: np.random.Generator):
        """Sample ``n`` droplets from the nozzle's fixed position.

        Returns ``(positions (n,3), velocities (n,3), radii (n,))``.
        """
        positions = np.tile(np.asarray(self.config.position, float), (n, 1))
        return self._emit_core(positions, rng)

    def emit_from(self, position, n: int, rng: np.random.Generator, *,
                  carriage_velocity=None):
        """Sample ``n`` droplets from an arbitrary ``position`` (m), optionally
        adding a nozzle ``carriage_velocity`` (m/s)."""
        positions = np.tile(np.asarray(position, float), (n, 1))
        cv = None if carriage_velocity is None else np.asarray(carriage_velocity, float)
        return self._emit_core(positions, rng, cv)

    def emit_path(self, moves, n: int, rng: np.random.Generator, *,
                  include_carriage_velocity: bool = True):
        """Sample ``n`` droplets spread along the spray-on (G1) segments of a path.

        Droplets are allocated to segments in proportion to their spray time
        (uniform deposition rate) and placed at uniformly random points along each
        segment. If ``include_carriage_velocity`` the segment's travel velocity is
        added to each droplet's launch velocity.
        """
        spray = [m for m in moves if m.spray_on and m.duration > 0.0 and m.length > 0.0]
        if n <= 0 or not spray:
            empty = np.zeros((0, 3))
            return empty, empty, np.zeros(0)

        durations = np.array([m.duration for m in spray])
        counts = _allocate_counts(n, durations / durations.sum())

        positions, carriage = [], []
        for move, count in zip(spray, counts):
            if count == 0:
                continue
            start = np.asarray(move.start, float)
            end = np.asarray(move.end, float)
            t = rng.uniform(0.0, 1.0, size=count)[:, None]
            positions.append(start + t * (end - start))
            seg_v = (end - start) / move.duration if include_carriage_velocity else np.zeros(3)
            carriage.append(np.tile(seg_v, (count, 1)))

        positions = np.vstack(positions)
        cv = np.vstack(carriage) if include_carriage_velocity else None
        return self._emit_core(positions, rng, cv)
