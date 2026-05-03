from __future__ import annotations

from .config import load_config
from .streamfunction import generate_streamfunction_member, save_member


def generate_ensemble(config_path: str) -> None:
    config = load_config(config_path)
    output_dir = config.paths.output_dir
    for member in range(config.stochastic_forcing.ensemble_size):
        result = generate_streamfunction_member(config, member)
        save_member(result, output_dir, config.paths.velocity_dir, member)
        print(f"Generated member {member + 1}/{config.stochastic_forcing.ensemble_size}")
