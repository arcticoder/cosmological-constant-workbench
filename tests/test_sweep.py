from __future__ import annotations

from ccw.mechanisms import CosmologyBackground
from ccw.sweep import run_sweep


def test_sweep_is_deterministic_and_sized() -> None:
    bg = CosmologyBackground()
    rows = run_sweep(
        mechanism="cpl",
        grid=[{"w0": -1.0, "wa": 0.0}, {"w0": -0.9, "wa": 0.0}],
        z_values=[0.0, 1.0, 2.0],
        bg=bg,
    )
    assert len(rows) == 2 * 3
    assert rows[0].mechanism in {"cpl_quintessence"}
