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
    rtol: float = 1e-6,
    atol: float = 1e-8,
) -> CoupledODEResult:
    """Solve coupled H(a) + scalar field equations.
    
    Parameters
    ----------
    bg : CosmologyBackground
        Background cosmology parameters.
    V_func : callable
        Potential V(φ) returning Energy Density [J/m³].
    V_prime_func : callable
        Derivative V'(φ) returning [J/m³].
    phi_0 : float
        Initial field value (dimensionless) at a=a_max.
    Pi_0 : float
        Initial momentum Π at a=a_max.
    a_min, a_max : float
        Integration range.
    n_eval, rtol, atol : 
        Solver parameters.
    """
    if a_min >= a_max:
        raise ValueError("a_min must be < a_max")
    
    # Physics Constants
    KAPPA = 8.0 * PI * G_M3_KG_S2  # SI units
    C2 = C_M_S**2
    
    # Background Setup
    h0_s_inv = h0_km_s_mpc_to_s_inv(bg.h0_km_s_mpc)
    rho_crit_0 = 3.0 * (h0_s_inv**2) / KAPPA  # Mass density [kg/m^3]
    
    rho_m_0 = bg.omega_m * rho_crit_0
    rho_r_0 = bg.omega_r * rho_crit_0
    
    # Time variable t = ln(a)
    t_min = math.log(a_min)
    t_max = math.log(a_max)
    
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        phi, Pi = y
        # Clamp Pi to avoid infinite loop in solver if it grows too large
        
        # Densities [kg/m^3]
        rho_m = rho_m_0 * math.exp(-3*t)
        rho_r = rho_r_0 * math.exp(-4*t)
        
        try:
            V_val = V_func(phi)
            V_prime_val = V_prime_func(phi)
        except (OverflowError, ValueError):
            return np.array([0.0, 0.0])
            
        rho_V = V_val / C2
        
        # Check Kinetic Dominance limit (Pi^2 < 6)
        if Pi**2 >= 5.99:
             # Prevent singularity by clamping effective Pi in denominator
             Pi_denom = 5.99
        else:
             Pi_denom = Pi**2
             
        # Friedmann Constraint
        # H^2 = (kappa/3) * (rho_fluid) / (1 - Pi^2/6)
        numerator = (KAPPA / 3.0) * (rho_m + rho_r + rho_V)
        denominator = 1.0 - Pi_denom / 6.0
        
        if numerator < 0: numerator = 1e-100
        if denominator <= 1e-4: denominator = 1e-4
             
        H_sq = numerator / denominator
             
        # EOMs
        # Term 1: Friction/Hubble Drag = -3 * Pi * (1 - Pi^2/6)
        term1 = -3.0 * Pi * (1.0 - Pi**2 / 6.0) 
        
        # Term 2: Fluid Coupling = (kappa / 2 H^2) * Pi * (rho_m + 4/3 rho_r)
        term2 = (KAPPA / (2.0 * max(H_sq, 1e-60))) * Pi * (rho_m + (4.0/3.0)*rho_r)
        
        # Term 3: Potential Gradient = - (kappa * V') / (c^2 * H^2)
        term3 = - (KAPPA * V_prime_val) / (C2 * max(H_sq, 1e-60))
        
        dphi_dt = Pi
        dPi_dt = term1 + term2 + term3
        
        return np.array([dphi_dt, dPi_dt])

    # Solve
    y0 = np.array([phi_0, Pi_0])
    sol = integrate.solve_ivp(
        rhs, 
        (t_max, t_min), 
        y0, 
        t_eval=np.linspace(t_max, t_min, n_eval),
        method='LSODA', 
        rtol=rtol, 
        atol=atol
    )
    
    # Process output
    t_out = sol.t[::-1]
    a_out = np.exp(t_out)
    phi_out = sol.y[0][::-1]
    Pi_out = sol.y[1][::-1]
    
    H_out = np.zeros_like(a_out)
    rho_DE_out = np.zeros_like(a_out)
    
    for i, (a, phi, Pi) in enumerate(zip(a_out, phi_out, Pi_out)):
        rho_m = rho_m_0 / (a**3)
        rho_r = rho_r_0 / (a**4)
        V_val = V_func(phi)
        rho_V = V_val / C2
        
        # Reconstruct H
        denom = 1.0 - min(Pi**2, 5.99) / 6.0
        numerator = (KAPPA / 3.0) * (rho_m + rho_r + rho_V)
        if numerator < 0: numerator = 0
        H_sq = max(1e-60, numerator / denom)
        H_val = math.sqrt(H_sq)
        
        H_out[i] = H_val
        
        # rho_DE is Energy Density [J/m^3]
        # rho_kin_mass = (3 H^2 / KAPPA) * (Pi^2 / 6) = H^2 Pi^2 / (2 KAPPA)
        rho_kin_mass = (H_sq * Pi**2) / (2.0 * KAPPA)
        rho_DE_mass = rho_kin_mass + rho_V
        rho_DE_out[i] = rho_DE_mass * C2
    
    return CoupledODEResult(
        a_grid=a_out,
        phi_grid=phi_out,
        Pi_grid=Pi_out,
        H_grid=H_out,
        rho_DE_grid=rho_DE_out,
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
    """Inverse power-law: V(φ) = M^(4+α) / φ^α."""
    if phi <= 0: return 1e100 # Soft barrier
    return M**(4 + alpha) / (phi**alpha)


def inverse_power_potential_prime(phi: float, M: float = 1e-3, alpha: float = 0.5) -> float:
    """Derivative of inverse power-law potential."""
    if phi <= 0: return -1e100 # Repulsive force
    return -alpha * M**(4 + alpha) / (phi**(alpha + 1))
