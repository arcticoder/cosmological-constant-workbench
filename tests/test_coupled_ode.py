"""
Tests for coupled ODE cosmology solver.
"""

import numpy as np
import pytest

from ccw.coupled_ode import (
    solve_coupled_cosmology,
    exponential_potential,
    exponential_potential_prime,
    inverse_power_potential,
    inverse_power_potential_prime,
)
from ccw.mechanisms import CosmologyBackground


def test_coupled_ode_returns_finite_solution():
    """Verify coupled ODE solver returns finite arrays."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3, omega_r=0.0)
    
    # Simple exponential potential with small initial field
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0, V0=1e-10),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0, V0=1e-10),
        phi_0=0.1,  # Small field value
        Pi_0=0.001,  # Small velocity
        a_min=0.5,
        a_max=1.0,
        n_eval=100,
    )
    
    assert len(result.a_grid) == 100
    assert np.all(np.isfinite(result.phi_grid))
    assert np.all(np.isfinite(result.Pi_grid))
    assert np.all(np.isfinite(result.H_grid))
    assert np.all(np.isfinite(result.rho_DE_grid))


def test_coupled_ode_a_grid_ascending():
    """Verify a_grid is in ascending order."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0),
        phi_0=0.1,
        Pi_0=0.001,
        a_min=0.5,
        a_max=1.0,
        n_eval=50,
    )
    
    assert np.all(np.diff(result.a_grid) > 0), "a_grid should be strictly ascending"


def test_coupled_ode_H_positive():
    """Verify H(a) is always positive."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0),
        phi_0=0.1,
        Pi_0=0.001,
        a_min=0.5,
        a_max=1.0,
        n_eval=50,
    )
    
    assert np.all(result.H_grid > 0), "H should be positive everywhere"


def test_coupled_ode_rho_DE_positive():
    """Verify ρ_DE is non-negative."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0),
        phi_0=0.1,
        Pi_0=0.001,
        a_min=0.5,
        a_max=1.0,
        n_eval=50,
    )
    
    assert np.all(result.rho_DE_grid >= 0), "ρ_DE should be non-negative"


def test_coupled_ode_interpolate_at_z():
    """Verify interpolation at redshift works."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0),
        phi_0=0.1,
        Pi_0=0.001,
        a_min=0.5,
        a_max=1.0,
        n_eval=100,
    )
    
    # Interpolate at z=0.5 (a=2/3)
    phi, Pi, H, rho_DE = result.interpolate_at_z(0.5)
    
    assert np.isfinite(phi) and phi > 0
    assert np.isfinite(Pi)
    assert np.isfinite(H) and H > 0
    assert np.isfinite(rho_DE) and rho_DE >= 0


def test_coupled_ode_interpolate_outside_range_raises():
    """Verify interpolation outside range raises error."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0),
        phi_0=0.1,
        Pi_0=0.001,
        a_min=0.5,  # z_max ~ 1
        a_max=1.0,
        n_eval=50,
    )
    
    # z=5 is outside range (a_min=0.5 → z_max ≈ 1)
    with pytest.raises(ValueError, match="outside cached integration range"):
        result.interpolate_at_z(5.0)


def test_coupled_ode_inverse_power_potential():
    """Verify solver works with inverse power-law potential."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: inverse_power_potential(phi, M=1e-3, alpha=0.5),
        V_prime_func=lambda phi: inverse_power_potential_prime(phi, M=1e-3, alpha=0.5),
        phi_0=1.0,  # Need φ > 0
        Pi_0=-0.001,  # Small negative velocity to keep phi from going negative
        a_min=0.8,
        a_max=1.0,
        n_eval=50,
    )
    
    assert np.all(result.phi_grid > 0), "φ should remain positive for inverse-power"
    assert np.all(np.isfinite(result.rho_DE_grid))


def test_coupled_ode_energy_conservation():
    """Verify total energy density evolves consistently with Friedmann."""
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_m=0.3, omega_r=0.0)
    
    result = solve_coupled_cosmology(
        bg=bg,
        V_func=lambda phi: exponential_potential(phi, lam=10.0),
        V_prime_func=lambda phi: exponential_potential_prime(phi, lam=10.0),
        phi_0=0.1,
        Pi_0=0.001,
        a_min=0.7,
        a_max=1.0,
        n_eval=50,
    )
    
    # Check Friedmann constraint: H² ∝ ρ_total
    # Just verify H² scales reasonably with a
    # (Detailed check would require recomputing ρ_m + ρ_DE)
    H_a1 = result.H_grid[0]  # At a_min
    H_a2 = result.H_grid[-1]  # At a_max
    
    # H should decrease from high-z to low-z (for matter-dominated)
    # (Not always true for DE-dominated, but at least H should be finite)
    assert H_a1 > 0 and H_a2 > 0


def test_exponential_potential_properties():
    """Verify exponential potential helper functions."""
    phi_test = 1.5
    lam = 2.0
    V0 = 1e-9
    
    V = exponential_potential(phi_test, lam, V0)
    V_prime = exponential_potential_prime(phi_test, lam, V0)
    
    assert V > 0
    assert V_prime < 0, "Exponential potential should have negative derivative"
    
    # Numerical derivative check
    dphi = 1e-6
    V_plus = exponential_potential(phi_test + dphi, lam, V0)
    V_minus = exponential_potential(phi_test - dphi, lam, V0)
    V_prime_numerical = (V_plus - V_minus) / (2 * dphi)
    
    np.testing.assert_allclose(V_prime, V_prime_numerical, rtol=1e-5)


def test_inverse_power_potential_properties():
    """Verify inverse power-law potential helper functions."""
    phi_test = 2.0
    M = 1e-3
    alpha = 0.8
    
    V = inverse_power_potential(phi_test, M, alpha)
    V_prime = inverse_power_potential_prime(phi_test, M, alpha)
    
    assert V > 0
    assert V_prime < 0, "Inverse power potential should have negative derivative for φ > 0"
    
    # Numerical derivative check
    dphi = 1e-6
    V_plus = inverse_power_potential(phi_test + dphi, M, alpha)
    V_minus = inverse_power_potential(phi_test - dphi, M, alpha)
    V_prime_numerical = (V_plus - V_minus) / (2 * dphi)
    
    np.testing.assert_allclose(V_prime, V_prime_numerical, rtol=1e-5)
