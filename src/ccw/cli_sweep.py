from __future__ import annotations

import argparse
import json
from pathlib import Path

from .mechanisms import CosmologyBackground
from .sweep import run_sweep, write_csv, write_json


def _parse_grid_json(grid_json: str) -> list[dict]:
    obj = json.loads(grid_json)
    if not isinstance(obj, list) or any(not isinstance(x, dict) for x in obj):
        raise ValueError("--grid-json must be a JSON array of objects")
    return obj


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ccw sweep runner (toy mechanisms)")
    p.add_argument("--mechanism", required=True, help="cpl | rvm | unimodular | sequestering")
    p.add_argument(
        "--z",
        default="0,1,2",
        help="Comma-separated list of redshifts (e.g. 0,0.5,1,2)",
    )
    p.add_argument(
        "--grid-json",
        default="[{}]",
        help=(
            "JSON array of param dicts, e.g. "
            "'[{\"w0\": -1.0, \"wa\": 0.2}, {\"w0\": -0.9, \"wa\": 0.0}]'"
        ),
    )
    p.add_argument("--h0", type=float, default=67.4)
    p.add_argument("--omega-lambda", type=float, default=0.6889)
    p.add_argument("--omega-m", type=float, default=0.3111)
    p.add_argument("--omega-r", type=float, default=0.0)
    p.add_argument("--omega-k", type=float, default=0.0)
    p.add_argument("--out-json", type=str, default="results/sweep.json")
    p.add_argument("--out-csv", type=str, default="results/sweep.csv")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    z_values = [float(s.strip()) for s in str(args.z).split(",") if s.strip()]
    grid = _parse_grid_json(args.grid_json)

    bg = CosmologyBackground(
        h0_km_s_mpc=args.h0,
        omega_lambda=args.omega_lambda,
        omega_m=args.omega_m,
        omega_r=args.omega_r,
        omega_k=args.omega_k,
    )

    rows = run_sweep(mechanism=args.mechanism, grid=grid, z_values=z_values, bg=bg)

    write_json(Path(args.out_json), rows)
    write_csv(Path(args.out_csv), rows)

    print(f"Wrote {len(rows)} rows to {args.out_json} and {args.out_csv}")
