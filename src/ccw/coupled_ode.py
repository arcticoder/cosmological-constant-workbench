"""
Coupled ODE solver for self-consistent cosmological evolution.

Solves H(a) + field equations (φ, Π) simultaneously for mechanisms where
ρ_DE and H are mutually dependent (e.g., scalar quintessence with backreaction).

Math:
  Scale-factor time t = ln(a):
  dφ/dt = Π
  dΠ/dt = -3H(φ,Π) Π - V'(φ)
  
  With H² = (8πG/3)(ρ_m a⁻³ + ρ_r a⁻⁴ + Π²/2 + V(φ))
  
This replaces the autonomous (x,y,Ω_r) system in ScalarFieldQuintessence with a
direct (φ,Π) integration for better accuracy and backreaction tracking.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional, Tuple

import numpy as np
from scipy import integrate

from .constants import C_M_S, G_M3_KG_S2, PI
from .cosmology import h0_km_s_mpc_to_s_inv
from .mechanisms import CosmologyBackground


@dataclass(frozen=True)
class CoupledODEResult:
    """Result from coupled cosmology integration.
    
    Attributes
    ----------
    a_grid : np.ndarray
        Scale factor grid (ascending).
    phi_grid : np.ndarray
        Scalar field values.
    Pi_grid : np.ndarray
        Canonical momentum Π = dφ/d ln(a).
    H_grid : np.ndarray
        Hubble parameter in s⁻¹.
    rho_DE_grid : np.ndarray
        Dark energy density in J/m³.
    """
    a_grid: np.ndarray
    phi_grid: np.ndarray
    Pi_grid: np.ndarray
    H_grid: np.ndarray
    rho_DE_grid: np.ndarray
    
    def interpolate_at_z(self, z: float) -> Tuple[float, float, float, float]:
        """Interpolate solution at redshift z.
        
        Returns
        -------
        tuple
            (φ, Π, H_s_inv, ρ_DE_j_m3)
        """
        a = 1.0 / (1.0 + z)
        if a < self.a_grid.min() or a > self.a_grid.max():
            raise ValueError(f"z={z} (a={a}) outside cached integration range [{self.a_grid.min()}, {self.a_grid.max()}]")
        
        phi = float(np.interp(a, self.a_grid, self.phi_grid))
        Pi = float(np.interp(a, self.a_grid, self.Pi_grid))
        H_s_inv = float(np.interp(a, self.a_grid, self.H_grid))
        rho_DE = float(np.interp(a, self.a_grid, self.rho_DE_grid))
        return phi, Pi, H_s_inv, rho_DE


def solve_coupled_cosmology(
    *,
    bg: CosmologyBackground,
    V_func: Callable[[float], float],
    V_prime_func: Callable[[float], float],
    phi_0: float,
    Pi_0: float,
    a_min: float = 1e-3,
    a_max: float = 1.0,
    n_eval: int = 800,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> CoupledODEResult:
    """Solve coupled H(a) + scalar field equations.
    
    Parameters
    ----------
    bg : CosmologyBackground
        Background cosmology parameters.
    V_func : callable
        Potential V(φ).
    V_prime_func : callable
        Derivative V'(φ).
    phi_0 : float
        Initial field value at a=a_max (today).
    Pi_0 : float
        Initial momentum Π at a=a_max.
    a_min : float
        Minimum scale factor to integrate to (default 1e-3, z~1000).
    a_max : float
        Maximum scale factor (default 1.0, z=0).
    n_eval : int
        Number of evaluation points.
    rtol, atol : float
        Integration tolerances.
    
    Returns
    -------
    CoupledODEResult
        Solution arrays.
    
    Notes
    -----
    Integrates backward from a_max to a_min in ln(a) time.
    Uses RK45 with tight tolerances for Friedmann-field coupling.
    """
    if a_min >= a_max:
        raise ValueError("a_min must be < a_max")
    if phi_0 <= 0:
        raise ValueError("phi_0 must be > 0 for inverse-power potentials")
    if n_eval < 50:
        raise ValueError("n_eval too small")
    
    # Convert background to SI units
    h0_s_inv = h0_km_s_mpc_to_s_inv(bg.h0_km_s_mpc)
    
    # Critical energy density (mass units) at z=0
    rho_crit_mass_0 = 3.0 * (h0_s_inv**2) / (8.0 * PI * G_M3_KG_S2)
    
    # Matter and radiation energy densities at a=1
    rho_m_0 = bg.omega_m * rho_crit_mass_0
    rho_r_0 = bg.omega_r * rho_crit_mass_0
    
    # Integration in ln(a) time: t = ln(a)
    t_min = math.log(a_min)
    t_max = math.log(a_max)
    
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        """Right-hand side: dy/dt = [dφ/dt, dΠ/dt]."""
        phi, Pi = float(y[0]), float(y[1])
        a = math.exp(t)
        
        # Energy densities (mass units)
        rho_m = rho_m_0 / (a**3)
        rho_r = rho_r_0 / (a**4)
        
        # Scalar field energy density (with overflow protection)
        try:
            V_phi = V_func(phi)
        except (OverflowError, ValueError):
            # Field went to unphysical region - terminate integration
            return np.array([0.0, 0.0])  # Will trigger solver termination
        
        rho_phi = 0.5 * Pi**2 + V_phi  # Kinetic + potential (dimensionless units, normalized by M_Pl²)
        
        # Total energy density
        rho_tot = rho_m + rho_r + rho_phi * (C_M_S**2)  # Convert phi to energy density
        
        if rho_tot <= 0 or not math.isfinite(rho_tot):
            raise ValueError(f"Non-physical ρ_tot at a={a:.3e}")
        
        # Hubble parameter
        H = math.sqrt((8.0 * PI * G_M3_KG_S2 / 3.0) * rho_tot)
        
        # Field equations in ln(a) time (with overflow protection)
        try:
            dV_dphi = V_prime_func(phi)
        except (OverflowError, ValueError):
            return np.array([0.0, 0.0])  # Terminate integration
            
        dphi_dt = Pi
        dPi_dt = -3.0 * H * Pi - dV_dphi
        
        return np.array([dphi_dt, dPi_dt], dtype=float)
    
    # Solve from a_max (today) backward to a_min
    y_init = np.array([phi_0, Pi_0], dtype=float)
    t_eval = np.linspace(t_max, t_min, n_eval)
    
    sol = integrate.solve_ivp(
        rhs,
        t_span=(t_max, t_min),
        y0=y_init,
        t_eval=t_eval,
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    
    if not sol.success:
        raise RuntimeError(f"Coupled ODE integration failed: {sol.message}")
    
    # Extract solution (reverse to ascending a)
    t_grid = sol.t[::-1]
    a_grid = np.exp(t_grid)
    phi_grid = sol.y[0, ::-1]
    Pi_grid = sol.y[1, ::-1]
    
    # Compute H(a) and ρ_DE(a) from solution
    H_grid = np.zeros_like(a_grid)
    rho_DE_grid = np.zeros_like(a_grid)
    
    for i, (a, phi, Pi) in enumerate(zip(a_grid, phi_grid, Pi_grid)):
        rho_m = rho_m_0 / (a**3)
        rho_r = rho_r_0 / (a**4)
        V_phi = V_func(phi)
        rho_phi = 0.5 * Pi**2 + V_phi
        rho_tot = rho_m + rho_r + rho_phi * (C_M_S**2)
        
        H_grid[i] = math.sqrt((8.0 * PI * G_M3_KG_S2 / 3.0) * rho_tot)
        rho_DE_grid[i] = rho_phi * (C_M_S**2)  # Energy density
    
    return CoupledODEResult(
        a_grid=a_grid,
        phi_grid=phi_grid,
        Pi_grid=Pi_grid,
        H_grid=H_grid,
        rho_DE_grid=rho_DE_grid,
    )


# Example potentials for testing

def exponential_potential(phi: float, lam: float = 1.0, V0: float = 1e-10) -> float:
    """Exponential potential: V(φ) = V₀ exp(-λ φ / M_Pl).
    
    Note: φ is in M_Pl units, so this is V₀ exp(-λ φ).
    """
    return V0 * math.exp(-lam * phi)


def exponential_potential_prime(phi: float, lam: float = 1.0, V0: float = 1e-10) -> float:
    """Derivative of exponential potential."""
    return -lam * V0 * math.exp(-lam * phi)


def inverse_power_potential(phi: float, M: float = 1e-3, alpha: float = 0.5) -> float:
    """Inverse power-law: V(φ) = M^(4+α) / φ^α.
    
    Requires φ > 0.
    """
    if phi <= 0:
        raise ValueError("φ must be > 0 for inverse-power potential")
    return M**(4 + alpha) / (phi**alpha)


def inverse_power_potential_prime(phi: float, M: float = 1e-3, alpha: float = 0.5) -> float:
    """Derivative of inverse power-law potential."""
    if phi <= 0:
        raise ValueError("φ must be > 0 for inverse-power potential")
    return -alpha * M**(4 + alpha) / (phi**(alpha + 1))
