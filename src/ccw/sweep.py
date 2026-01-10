from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .baseline import BaselineInputs, compute_baseline
from .mechanisms import CPLQuintessence, CosmologyBackground, RunningVacuumRVM, SequesteringToy, UnimodularBookkeeping


@dataclass(frozen=True)
class SweepRow:
    mechanism: str
    params: Dict[str, Any]
    z: float
    rho_de_j_m3: float
    w_de: Optional[float]


def evaluate_mechanism(name: str, params: Dict[str, Any], z: float, bg: CosmologyBackground) -> SweepRow:
    key = name.strip().lower()

    if key in {"cpl", "cpl_quintessence"}:
        mech = CPLQuintessence(w0=float(params.get("w0", -1.0)), wa=float(params.get("wa", 0.0)))
    elif key in {"rvm", "running_vacuum", "running_vacuum_rvm"}:
        mech = RunningVacuumRVM(nu=float(params.get("nu", 0.0)))
    elif key in {"unimodular", "unimodular_bookkeeping"}:
        mech = UnimodularBookkeeping()
    elif key in {"sequestering", "sequestering_toy"}:
        mech = SequesteringToy(delta_rho_j_m3=float(params.get("delta_rho_j_m3", 0.0)))
    else:
        raise ValueError(f"Unknown mechanism: {name}")

    out = mech.evaluate(z=z, bg=bg).result
    return SweepRow(mechanism=mech.name, params=params, z=z, rho_de_j_m3=out.rho_de_j_m3, w_de=out.w_de)


def run_sweep(
    *,
    mechanism: str,
    grid: Iterable[Dict[str, Any]],
    z_values: Iterable[float],
    bg: CosmologyBackground,
) -> List[SweepRow]:
    rows: List[SweepRow] = []
    for params in grid:
        for z in z_values:
            rows.append(evaluate_mechanism(mechanism, params, float(z), bg))
    return rows


def write_json(path: Path, rows: List[SweepRow]) -> None:
    payload = [asdict(r) for r in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def write_csv(path: Path, rows: List[SweepRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mechanism", "z", "rho_de_j_m3", "w_de", "params_json"])
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "mechanism": r.mechanism,
                    "z": r.z,
                    "rho_de_j_m3": r.rho_de_j_m3,
                    "w_de": r.w_de,
                    "params_json": json.dumps(r.params, sort_keys=True),
                }
            )
