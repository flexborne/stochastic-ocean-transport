#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.circulation.generate_ensemble import generate_ensemble


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate stochastic Stommel stream-function ensemble.")
    parser.add_argument("--config", default="configs/demo.yaml", help="Path to YAML configuration file.")
    args = parser.parse_args()
    generate_ensemble(args.config)


if __name__ == "__main__":
    main()
