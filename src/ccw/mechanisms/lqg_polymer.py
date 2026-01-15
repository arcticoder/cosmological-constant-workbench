"""LQG polymer (holonomy) cosmology: effective dark-energy bookkeeping.

This module provides a lightweight, algebraic way to explore loop-quantum-cosmology
(LQC/LQG) inspired holonomy corrections without relying on the blocked coupled ODE
solver (J.22).

Core idea (standard LQC form):

  H(z)^2 = (8πG/3c^2) ρ(z) [1 - ρ(z)/ρ_c]

where ρ(z) is the total energy density (J/m^3) and ρ_c is a critical density
(typically ~0.41 ρ_Planck in the "improved dynamics" scheme).

To integrate with the workbench's existing interface (mechanisms provide ρ_DE(z)),
we rewrite this as an *effective* dark-energy density:

  ρ_eff(z) = ρ(z) + ρ_DE,eff(z)
  ρ_DE,eff(z) := - ρ(z)^2 / ρ_c

so that the usual Friedmann equation recovers the LQC correction when using
ρ_total = ρ_m + ρ_r + ρ_DE,eff.

Important:
- This effective ρ_DE,eff is generally NEGATIVE and highly redshift-dependent.
- As implemented here, it is a diagnostic / exploratory mechanism. It is not
  claimed to solve the cosmological constant problem.

We include a toy parameter `mu0_factor` that rescales ρ_c via ρ_c ∝ 1/μ0^2.
In physical LQC, μ0 is set by the area gap ~ l_Pl; here we allow scans as an
exploratory knob to see what scales would be required to impact late-time
cosmology.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .base import CosmologyBackground, MechanismOutput, MechanismResult, ensure_z_nonnegative
from ..constants import C_M_S, G_M3_KG_S2, HBAR_J_S, PI


def planck_energy_density_j_m3() -> float:
    """Planck energy density in J/m^3.

    Using the standard combination:
      ρ_Pl,energy = c^7 / (ħ G^2)

    Units check: J/m^3.
    """
    return (C_M_S**7) / (HBAR_J_S * (G_M3_KG_S2**2))


def lqc_critical_density_j_m3(
    *,
    rho_c_over_rho_pl: float = 0.41,
    mu0_factor: float = 1.0,
) -> float:
    """Toy critical density ρ_c in J/m^3.

    Parameters
    ----------
    rho_c_over_rho_pl:
        Baseline ratio ρ_c / ρ_Pl. In common LQC conventions this is ~0.41.
    mu0_factor:
        Dimensionless scale factor for μ0 relative to a Planck-scale choice.

        In LQC, ρ_c ∝ 1/μ0^2, so here we model
          ρ_c(mu0_factor) = ρ_c(1) / mu0_factor^2.

    Notes
    -----
    This is intentionally a *scan parameter*. Large mu0_factor values correspond
    to lowering the bounce density dramatically, which is not standard.
    """
    if rho_c_over_rho_pl <= 0:
        raise ValueError("rho_c_over_rho_pl must be positive")
    if mu0_factor <= 0:
        raise ValueError("mu0_factor must be positive")

    rho_pl = planck_energy_density_j_m3()
    return (rho_c_over_rho_pl * rho_pl) / (mu0_factor**2)


@dataclass(frozen=True)
class LQGPolymerCosmology:
    """Effective dark-energy mechanism capturing LQC holonomy corrections.

    Parameters
    ----------
    rho_c_over_rho_pl:
        Baseline critical density ratio. Default ~0.41.
    mu0_factor:
        Rescales ρ_c via ρ_c ∝ 1/μ0^2.
    include_lambda:
        If True, adds the observed ρ_Λ,0 as a separate constant component.
        This is useful for testing how polymer corrections would *perturb* ΛCDM,
        but it is not a derivation of Λ.
    """

    rho_c_over_rho_pl: float = 0.41
    mu0_factor: float = 1.0
    include_lambda: bool = False

    name: str = "lqg_polymer"

    def __post_init__(self) -> None:
        if self.rho_c_over_rho_pl <= 0:
            raise ValueError("rho_c_over_rho_pl must be positive")
        if self.mu0_factor <= 0:
            raise ValueError("mu0_factor must be positive")

    def describe_assumptions(self) -> str:
        base = (
            "LQC/LQG holonomy correction treated as an effective density "
            "ρ_DE,eff(z) = -ρ(z)^2/ρ_c with ρ_c scaled as ρ_c ∝ 1/μ0^2. "
            "This is exploratory bookkeeping; ρ_DE,eff is typically negative."
        )
        if self.include_lambda:
            base += " Includes an explicit constant ρ_Λ,0 added separately."
        return base

    def rho_c_j_m3(self) -> float:
        return lqc_critical_density_j_m3(rho_c_over_rho_pl=self.rho_c_over_rho_pl, mu0_factor=self.mu0_factor)

    def _rho_crit0_j_m3(self, bg: CosmologyBackground) -> float:
        if bg.omega_lambda <= 0:
            raise ValueError("omega_lambda must be > 0 to infer rho_crit0 from rho_lambda0")
        return bg.rho_lambda0_j_m3 / bg.omega_lambda

    def _rho_m_j_m3(self, z: float, bg: CosmologyBackground) -> float:
        rho_crit0 = self._rho_crit0_j_m3(bg)
        return bg.omega_m * rho_crit0 * (1.0 + z) ** 3

    def _rho_r_j_m3(self, z: float, bg: CosmologyBackground) -> float:
        rho_crit0 = self._rho_crit0_j_m3(bg)
        return bg.omega_r * rho_crit0 * (1.0 + z) ** 4

    def effective_rho_de_j_m3(self, z: float, bg: CosmologyBackground) -> float:
        """Return ρ_DE,eff(z) in J/m^3.

        Uses ρ(z) = ρ_m(z) + ρ_r(z) for the holonomy correction.
        """
        ensure_z_nonnegative(z)
        rho_m = self._rho_m_j_m3(z, bg)
        rho_r = self._rho_r_j_m3(z, bg)
        rho = rho_m + rho_r

        rho_c = self.rho_c_j_m3()
        rho_de_eff = -(rho**2) / rho_c

        if self.include_lambda:
            rho_de_eff += bg.rho_lambda0_j_m3

        return float(rho_de_eff)

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:
        ensure_z_nonnegative(z)
        rho_de = self.effective_rho_de_j_m3(z, bg)

        return MechanismOutput(
            result=MechanismResult(z=z, rho_de_j_m3=rho_de, w_de=None),
            assumptions=self.describe_assumptions(),
        )


def lambda_eff_m_minus2_from_h0_omega(h0_km_s_mpc: float, omega_lambda: float) -> float:
    """Convenience: Λ_eff for a flat FRW with (H0, ΩΛ).

    Λ = 3 ΩΛ H0^2 / c^2.

    Returns Λ in m^-2.
    """
    if h0_km_s_mpc <= 0:
        raise ValueError("h0_km_s_mpc must be positive")
    if omega_lambda < 0:
        raise ValueError("omega_lambda must be >= 0")

    # H0: km/s/Mpc -> s^-1
    # 1 Mpc = 3.085677581e22 m
    mpc_m = 3.085677581e22
    h0_s_inv = (h0_km_s_mpc * 1e3) / mpc_m

    return float(3.0 * omega_lambda * (h0_s_inv**2) / (C_M_S**2))
