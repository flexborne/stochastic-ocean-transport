from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class DomainConfig:
    length_x_m: float
    length_y_m: float
    grid_x: int
    grid_y: int


@dataclass(frozen=True)
class StommelConfig:
    depth_m: float
    friction_coefficient_m_s: float
    beta_m_inv_s_inv: float
    wind_stress_curl_amplitude_m2_s2: float
    fourier_modes: int


@dataclass(frozen=True)
class StochasticForcingConfig:
    enabled: bool
    relative_std: float
    correlation_radius_m: float
    ensemble_size: int
    random_seed: int | None = None


@dataclass(frozen=True)
class PathConfig:
    input_dir: Path
    velocity_dir: Path
    output_dir: Path
    reference_dir: Path


@dataclass(frozen=True)
class AppConfig:
    domain: DomainConfig
    stommel: StommelConfig
    stochastic_forcing: StochasticForcingConfig
    paths: PathConfig
    raw: dict[str, Any]


def load_config(path: str | Path) -> AppConfig:
    cfg_path = Path(path)
    with cfg_path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)

    base = cfg_path.parent.parent
    paths_raw = raw["paths"]
    paths = PathConfig(
        input_dir=(base / paths_raw["input_dir"]).resolve(),
        velocity_dir=(base / paths_raw["velocity_dir"]).resolve(),
        output_dir=(base / paths_raw["output_dir"]).resolve(),
        reference_dir=(base / paths_raw["reference_dir"]).resolve(),
    )
    return AppConfig(
        domain=DomainConfig(**raw["domain"]),
        stommel=StommelConfig(**raw["stommel"]),
        stochastic_forcing=StochasticForcingConfig(**raw["stochastic_forcing"]),
        paths=paths,
        raw=raw,
    )
