from __future__ import annotations

import numpy as np
from PyAstronomy import pyasl


def exponentially_correlated_noise(
    n_points: int,
    correlation_radius_m: float,
    target_std: float,
    seed: int | None = None,
) -> np.ndarray:
    if n_points <= 1:
        raise ValueError("n_points must be greater than 1")
    if correlation_radius_m <= 0.0:
        raise ValueError("correlation_radius_m must be positive")

    if seed is not None:
        np.random.seed(seed)

    deviation, _, _ = pyasl.expCorrRN(
        n_points,
        correlation_radius_m,
        mean=0.0,
        std=1.0,
        fullOut=True,
    )
    deviation = np.asarray(deviation, dtype=float)

    std = float(np.std(deviation))
    if std == 0.0:
        raise ValueError("Generated noise has zero standard deviation.")
    return deviation * (target_std / std)
