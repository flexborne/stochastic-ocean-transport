from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import integrate

from .config import AppConfig
from .stochastic_forcing import exponentially_correlated_noise


@dataclass(frozen=True)
class StreamfunctionResult:
    psi: np.ndarray
    u: np.ndarray
    v: np.ndarray
    deviation: np.ndarray
    fourier_coefficients: np.ndarray
    gamma: np.ndarray


def _lambda(k: int, length_y: float) -> float:
    return math.pi * k / length_y


def _y_mode(y: np.ndarray | float, k: int, length_y: float) -> np.ndarray | float:
    return np.sin(_lambda(k, length_y) * y)


def _compute_velocity_from_streamfunction(psi: np.ndarray, length_x: float, length_y: float) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = psi.shape
    dy = length_y / max(ny - 1, 1)
    dx = length_x / max(nx - 1, 1)
    dpsi_dy, dpsi_dx = np.gradient(psi, dy, dx, edge_order=2)
    return dpsi_dy, -dpsi_dx


def generate_streamfunction_member(config: AppConfig, member_index: int) -> StreamfunctionResult:
    domain = config.domain
    stommel = config.stommel
    forcing = config.stochastic_forcing

    length_x = domain.length_x_m
    length_y = domain.length_y_m
    modes = stommel.fourier_modes
    forcing_amplitude = stommel.wind_stress_curl_amplitude_m2_s2

    len_d = math.ceil(length_y)
    target_std = forcing.relative_std * forcing_amplitude
    seed = None if forcing.random_seed is None else forcing.random_seed + member_index
    deviation = exponentially_correlated_noise(
        len_d,
        forcing.correlation_radius_m,
        target_std,
        seed=seed,
    )

    def wind_stress(y: float) -> float:
        idx = math.floor(y * (len_d - 1) / length_y)
        idx = max(0, min(len_d - 1, idx))
        return (forcing_amplitude + deviation[idx]) * math.cos(math.pi * y / length_y)

    coeffs: list[float] = []
    gamma: list[float] = []
    for k in range(1, modes + 1):
        res = lambda s, kk=k: wind_stress(s) * math.cos(_lambda(kk, length_y) * s)
        value = integrate.quad(res, 0.0, length_y, limit=100)[0] * 2.0 / length_y
        coeffs.append(value)
        gamma.append(value * math.pi * k / (stommel.friction_coefficient_m_s * length_y))
    coeffs_np = np.asarray(coeffs, dtype=float)
    gamma_np = np.asarray(gamma, dtype=float)

    alpha = 0.0

    def x_component(x: float, k: int) -> float:
        lam = _lambda(k, length_y)
        a_exp = -alpha * (stommel.depth_m / stommel.friction_coefficient_m_s) / 2.0 + math.sqrt(
            (alpha * (stommel.depth_m / stommel.friction_coefficient_m_s)) ** 2 / 4.0 + lam * lam
        )
        b_exp = -alpha * (stommel.depth_m / stommel.friction_coefficient_m_s) / 2.0 - math.sqrt(
            (alpha * (stommel.depth_m / stommel.friction_coefficient_m_s)) ** 2 / 4.0 + lam * lam
        )
        n_value = gamma_np[k - 1] * (length_y / (math.pi * k)) ** 2
        denominator = math.exp(a_exp * length_x) - math.exp(b_exp * length_x)
        p = n_value * (1.0 - math.exp(b_exp * length_x)) / denominator
        q = n_value - p
        return p * math.exp(a_exp * x) + q * math.exp(b_exp * x)

    x_values = np.linspace(0.0, length_x, domain.grid_x)
    y_values = np.linspace(0.0, length_y, domain.grid_y)
    psi = np.zeros((domain.grid_y, domain.grid_x), dtype=float)
    for iy, y in enumerate(y_values):
        for ix, x in enumerate(x_values):
            homogeneous = sum(_y_mode(y, k, length_y) * x_component(x, k) for k in range(1, modes + 1))
            particular = sum(
                gamma_np[k - 1] * (length_y / (math.pi * k)) ** 2 * _y_mode(y, k, length_y)
                for k in range(1, modes + 1)
            )
            psi[iy, ix] = homogeneous - particular

    u, v = _compute_velocity_from_streamfunction(psi, length_x, length_y)
    return StreamfunctionResult(psi=psi, u=u, v=v, deviation=deviation, fourier_coefficients=coeffs_np, gamma=gamma_np)


def save_member(result: StreamfunctionResult, output_dir: Path, velocity_dir: Path, member_index: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    velocity_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_dir / f"grid{member_index}.txt", result.psi, fmt="%.8e")
    np.savetxt(output_dir / f"deviation{member_index}.txt", result.deviation.reshape(1, -1), fmt="%.8e")
    np.savetxt(velocity_dir / f"{member_index}_u.txt", result.u, fmt="%.8e")
    np.savetxt(velocity_dir / f"{member_index}_v.txt", result.v, fmt="%.8e")
