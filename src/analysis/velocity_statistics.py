from __future__ import annotations

from pathlib import Path

import numpy as np

from .statistics import correlation_with_reference, ensemble_dispersion, ensemble_mean, read_grid_series
from src.circulation.config import load_config


def spectral_derivative(field: np.ndarray, spacing: float, axis: int) -> np.ndarray:
    n = field.shape[axis]
    kappa = 2.0 * np.pi * np.fft.fftfreq(n, d=spacing)
    shape = [1] * field.ndim
    shape[axis] = n
    derivative_hat = 1j * kappa.reshape(shape) * np.fft.fft(field, axis=axis)
    return np.real(np.fft.ifft(derivative_hat, axis=axis))


def compute_velocity_from_streamfunction(psi: np.ndarray, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray]:
    # Oceanographic convention used in the paper: u = dpsi/dy, v = -dpsi/dx.
    u = spectral_derivative(psi, dy, axis=1)
    v = -spectral_derivative(psi, dx, axis=2)
    return u, v


def analyze_velocity(config_path: str) -> None:
    config = load_config(config_path)
    count = config.stochastic_forcing.ensemble_size
    fields = read_grid_series(config.paths.output_dir, "grid", count)
    dx = config.domain.length_x_m / max(config.domain.grid_x - 1, 1)
    dy = config.domain.length_y_m / max(config.domain.grid_y - 1, 1)
    u, v = compute_velocity_from_streamfunction(fields, dx, dy)

    output = config.paths.output_dir
    output.mkdir(parents=True, exist_ok=True)
    for name, values in {
        "u_mean.txt": ensemble_mean(u),
        "v_mean.txt": ensemble_mean(v),
        "u_dispersion.txt": ensemble_dispersion(u),
        "v_dispersion.txt": ensemble_dispersion(v),
        "u_correlation_zero.txt": correlation_with_reference(u, (1, 1)),
        "v_correlation_zero.txt": correlation_with_reference(v, (1, 1)),
        "u_correlation_center.txt": correlation_with_reference(u, (u.shape[1] // 2, u.shape[2] // 2)),
        "v_correlation_center.txt": correlation_with_reference(v, (v.shape[1] // 2, v.shape[2] // 2)),
    }.items():
        np.savetxt(output / name, values, fmt="%.5f")
