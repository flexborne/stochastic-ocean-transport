#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.analysis.velocity_statistics import analyze_velocity


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute velocity and ensemble statistics from generated stream functions.")
    parser.add_argument("--config", default="configs/demo.yaml", help="Path to YAML configuration file.")
    args = parser.parse_args()
    analyze_velocity(args.config)


if __name__ == "__main__":
    main()
