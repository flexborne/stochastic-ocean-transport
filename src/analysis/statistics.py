from __future__ import annotations

from pathlib import Path

import numpy as np


def read_grid_series(directory: Path, prefix: str, count: int) -> np.ndarray:
    fields = []
    for member in range(count):
        path = directory / f"{prefix}{member}.txt"
        if not path.exists():
            raise FileNotFoundError(f"Missing ensemble file: {path}")
        fields.append(np.loadtxt(path))
    return np.asarray(fields, dtype=float)


def ensemble_mean(fields: np.ndarray) -> np.ndarray:
    return np.mean(fields, axis=0)


def ensemble_dispersion(fields: np.ndarray) -> np.ndarray:
    mean = ensemble_mean(fields)
    return np.sum((fields - mean) ** 2, axis=0)


def correlation_with_reference(fields: np.ndarray, reference_index: tuple[int, int]) -> np.ndarray:
    mean = ensemble_mean(fields)
    dispersion = ensemble_dispersion(fields)
    iy, ix = reference_index
    ref_series = fields[:, iy, ix]
    ref_mean = mean[iy, ix]
    ref_dispersion = dispersion[iy, ix]
    result = np.zeros_like(mean)
    for y in range(fields.shape[1]):
        for x in range(fields.shape[2]):
            denom = np.sqrt(dispersion[y, x] * ref_dispersion)
            if denom != 0.0:
                result[y, x] = np.sum((fields[:, y, x] - mean[y, x]) * (ref_series - ref_mean)) / denom
    return result
