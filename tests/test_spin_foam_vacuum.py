"""Tests for toy spin-foam vacuum mechanism (Phase K.27 scaffolding)."""

import math

import pytest

from ccw.mechanisms import CosmologyBackground, SpinFoamVacuum
from ccw.mechanisms.spin_foam_vacuum import hubble_entropy_log10, planck_length_m


def test_planck_length_reasonable_magnitude():
    lp = planck_length_m()
    assert math.isfinite(lp)
    assert 1e-36 < lp < 1e-34


def test_hubble_entropy_log10_reasonable_for_planck_2018_h0():
    # Order-of-magnitude sanity check: log10(S_H) ~ 122-123 today.
    log10_s = hubble_entropy_log10(67.4)
    assert math.isfinite(log10_s)
    assert 121.0 < log10_s < 124.0


def test_spin_foam_vacuum_scaling_and_equation_of_state():
    bg = CosmologyBackground()

    mech_122 = SpinFoamVacuum(avg_amplitude=1.0, suppression_exponent_log10=122.0)
    mech_121 = SpinFoamVacuum(avg_amplitude=1.0, suppression_exponent_log10=121.0)

    rho_122 = mech_122.evaluate(z=0.0, bg=bg).result.rho_de_j_m3
    rho_121 = mech_121.evaluate(z=0.0, bg=bg).result.rho_de_j_m3

    assert rho_122 > 0
    assert rho_121 > rho_122

    # Reducing the suppression exponent by 1 should increase rho by ~10x.
    ratio = rho_121 / rho_122
    assert 9.5 < ratio < 10.5

    out = mech_122.evaluate(z=0.0, bg=bg)
    assert out.result.w_de == -1.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"gamma": 0.0},
        {"gamma": -1.0},
        {"avg_amplitude": 0.0},
        {"avg_amplitude": -1.0},
        {"suppression_exponent_log10": -0.1},
    ],
)
def test_spin_foam_vacuum_rejects_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        SpinFoamVacuum(**kwargs)
