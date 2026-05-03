from __future__ import annotations

import numpy as np


def exponentially_correlated_noise(
    n_points: int,
    correlation_radius_m: float,
    target_std: float,
    seed: int | None = None,
) -> np.ndarray:
    """Generate a zero-mean exponentially correlated random process.

    The process is generated as an AR(1) sequence whose correlation decays
    approximately as exp(-r / correlation_radius_m) on a unit-spaced grid.
    The resulting series is centered and rescaled to the requested standard
    deviation.
    """
    if n_points <= 1:
        raise ValueError("n_points must be greater than 1")
    if correlation_radius_m <= 0.0:
        raise ValueError("correlation_radius_m must be positive")

    rng = np.random.default_rng(seed)
    phi = float(np.exp(-1.0 / correlation_radius_m))
    innovation_std = float(np.sqrt(max(0.0, 1.0 - phi * phi)))
    values = np.empty(n_points, dtype=float)
    values[0] = rng.normal()
    for idx in range(1, n_points):
        values[idx] = phi * values[idx - 1] + innovation_std * rng.normal()

    values -= float(np.mean(values))
    std = float(np.std(values))
    if std == 0.0:
        raise ValueError("Generated noise has zero standard deviation.")
    return values * (target_std / std)
