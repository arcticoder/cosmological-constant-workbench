from __future__ import annotations

import math

from ccw.mechanisms import CPLQuintessence, CosmologyBackground, RunningVacuumRVM, SequesteringToy, UnimodularBookkeeping


def test_cpl_reduces_to_lambda_when_w0_minus1_wa0() -> None:
    bg = CosmologyBackground()
    m = CPLQuintessence(w0=-1.0, wa=0.0)

    r0 = m.evaluate(0.0, bg).result.rho_de_j_m3
    r2 = m.evaluate(2.0, bg).result.rho_de_j_m3

    assert math.isfinite(r0)
    assert math.isfinite(r2)
    assert abs(r0 - bg.rho_lambda0_j_m3) / bg.rho_lambda0_j_m3 < 1e-12
    assert abs(r2 - bg.rho_lambda0_j_m3) / bg.rho_lambda0_j_m3 < 1e-12


def test_unimodular_stub_is_constant() -> None:
    bg = CosmologyBackground()
    m = UnimodularBookkeeping()
    assert m.evaluate(0.0, bg).result.rho_de_j_m3 == bg.rho_lambda0_j_m3
    assert m.evaluate(5.0, bg).result.rho_de_j_m3 == bg.rho_lambda0_j_m3


def test_sequestering_delta_rho_shifts_constant() -> None:
    bg = CosmologyBackground()
    m = SequesteringToy(delta_rho_j_m3=1e-12)
    out = m.evaluate(1.0, bg).result
    assert abs(out.rho_de_j_m3 - (bg.rho_lambda0_j_m3 + 1e-12)) < 1e-18


def test_running_vacuum_nu0_is_reasonable_at_z0() -> None:
    bg = CosmologyBackground()
    m = RunningVacuumRVM(nu=0.0)
    out0 = m.evaluate(0.0, bg).result

    assert math.isfinite(out0.rho_de_j_m3)
    assert out0.rho_de_j_m3 > 0
