"""Droplet emitter: samples initial positions, velocities and radii."""

from __future__ import annotations

import numpy as np

from .config import NozzleConfig


def _orthonormal_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u, v, w) with w == normalised axis and u, v spanning its plane."""
    w = axis / np.linalg.norm(axis)
    # Pick a helper vector that is not parallel to w.
    helper = np.array([1.0, 0.0, 0.0]) if abs(w[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, w)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)
    return u, v, w


class Nozzle:
    """Emits droplets into a cone around a spray axis."""

    def __init__(self, config: NozzleConfig):
        self.config = config
        self._u, self._v, self._w = _orthonormal_basis(np.asarray(config.direction, float))

    def emit(self, n: int, rng: np.random.Generator):
        """Sample ``n`` droplets.

        Returns
        -------
        positions : (n, 3) array of launch positions (m)
        velocities : (n, 3) array of launch velocities (m/s)
        radii : (n,) array of droplet radii (m)
        """
        cfg = self.config

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

        # Speeds: positive, normally distributed around the exit speed.
        speeds = rng.normal(cfg.exit_speed, cfg.exit_speed * cfg.speed_spread, size=n)
        speeds = np.clip(speeds, 0.0, None)
        velocities = dirs * speeds[:, None]

        # Radii: log-normal so the distribution is positive and right-skewed.
        sigma = cfg.radius_spread
        mu = np.log(cfg.mean_radius) - 0.5 * sigma**2  # keeps the mean at mean_radius
        radii = rng.lognormal(mean=mu, sigma=sigma, size=n)

        positions = np.tile(np.asarray(cfg.position, float), (n, 1))
        return positions, velocities, radii
