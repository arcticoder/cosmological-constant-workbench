import math

from ccw.mechanisms import CosmologyBackground
from ccw.mechanisms.lqg_polymer import (
    LQGPolymerCosmology,
    lqc_critical_density_j_m3,
    planck_energy_density_j_m3,
)


def test_planck_energy_density_positive_and_finite():
    rho_pl = planck_energy_density_j_m3()
    assert math.isfinite(rho_pl)
    assert rho_pl > 0


def test_lqc_critical_density_scales_with_mu0_factor():
    rho_c_1 = lqc_critical_density_j_m3(mu0_factor=1.0)
    rho_c_10 = lqc_critical_density_j_m3(mu0_factor=10.0)
    # rho_c ∝ 1/mu0^2
    assert math.isclose(rho_c_10, rho_c_1 / 100.0, rel_tol=1e-12, abs_tol=0.0)


def test_polymer_correction_negligible_for_planck_scale_mu0():
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)
    mech = LQGPolymerCosmology(mu0_factor=1.0)

    rho_de0 = mech.effective_rho_de_j_m3(0.0, bg)
    # Expect an utterly negligible correction at late times for ρ_c ~ ρ_Pl.
    ratio = abs(rho_de0) / bg.rho_lambda0_j_m3
    assert ratio < 1e-80


def test_polymer_effect_scales_like_mu0_squared_at_fixed_z():
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)
    z = 0.0

    mech_a = LQGPolymerCosmology(mu0_factor=1e30)
    mech_b = LQGPolymerCosmology(mu0_factor=1e31)

    rho_a = abs(mech_a.effective_rho_de_j_m3(z, bg))
    rho_b = abs(mech_b.effective_rho_de_j_m3(z, bg))

    # Since rho_de ~ -rho^2 / rho_c and rho_c ∝ 1/mu0^2, we expect rho_de ∝ mu0^2.
    # So stepping mu0 by 10x should scale the magnitude by ~100x.
    assert math.isclose(rho_b / rho_a, 100.0, rel_tol=0.05, abs_tol=0.0)


def test_mechanism_evaluate_returns_output_with_optional_w():
    bg = CosmologyBackground()
    mech = LQGPolymerCosmology(mu0_factor=1.0)

    out = mech.evaluate(z=1.0, bg=bg)
    assert out.result.z == 1.0
    assert math.isfinite(out.result.rho_de_j_m3)
    assert out.result.w_de is None
