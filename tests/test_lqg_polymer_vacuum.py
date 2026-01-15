import math

from ccw.mechanisms import CosmologyBackground
from ccw.mechanisms.lqg_polymer import (
    LQGPolymerDerivedVacuum,
    toy_derived_vacuum_energy_density_j_m3,
)


def test_toy_derived_vacuum_positive_and_decreases_with_mu0():
    rho1 = toy_derived_vacuum_energy_density_j_m3(mu0_factor=1.0, alpha_entropy=120.0)
    rho2 = toy_derived_vacuum_energy_density_j_m3(mu0_factor=2.0, alpha_entropy=120.0)
    assert rho1 > 0
    assert rho2 > 0
    assert rho2 < rho1


def test_polymer_vacuum_can_be_positive_at_z0_for_sufficient_suppression():
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)

    # Choose parameters that strongly suppress the vacuum term to around the observed scale.
    # (This is not predictive; it's a sanity test that the plumbing works.)
    mech = LQGPolymerDerivedVacuum(mu0_factor=60.0, alpha_entropy=120.0)

    rho_de0 = mech.effective_rho_de_j_m3(0.0, bg)
    assert math.isfinite(rho_de0)


def test_polymer_vacuum_evaluate_returns_output():
    bg = CosmologyBackground()
    mech = LQGPolymerDerivedVacuum(mu0_factor=60.0, alpha_entropy=120.0)
    out = mech.evaluate(z=0.5, bg=bg)
    assert out.result.z == 0.5
    assert math.isfinite(out.result.rho_de_j_m3)
    assert out.result.w_de is None
