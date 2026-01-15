"""Toy spin-foam vacuum mechanism (Phase K.27 scaffolding).

This is NOT an EPRL/BC vertex-amplitude computation. It is a controlled toy
parameterization intended to answer a narrow question:

Can a horizon-area-related suppression of order ~10^-122 produce the observed
late-time vacuum energy scale without introducing a *new* independent hierarchy?

Model (toy):

  ρ_Λ,eff = ρ_Pl * Ā * 10^{-X}

where:
- ρ_Pl is the Planck energy density in J/m^3
- Ā is a dimensionless coarse-grained amplitude factor (expected O(1) in the
  optimistic reading; treated as a scan knob here)
- X is a dimensionless suppression exponent (base-10), motivated by the fact
  that the Hubble-horizon entropy today is ~10^122 in Planck units.

We provide helpers to compute the *suggested* entropy scale from H0:

  S_H = A_H / (4 l_P^2),  A_H = 4π (c/H0)^2

and return log10(S_H), which is typically ~122–123.

The mechanism returns constant w=-1 (pure vacuum energy bookkeeping).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .base import CosmologyBackground, MechanismOutput, MechanismResult, ensure_z_nonnegative
from ..constants import C_M_S, G_M3_KG_S2, HBAR_J_S, PI
from ..cosmology import h0_km_s_mpc_to_s_inv
from .lqg_polymer import planck_energy_density_j_m3


def planck_length_m() -> float:
    """Planck length ℓ_P in meters."""
    return math.sqrt((HBAR_J_S * G_M3_KG_S2) / (C_M_S**3))


def hubble_horizon_area_m2(h0_km_s_mpc: float) -> float:
    """Hubble-horizon area A_H = 4π (c/H0)^2 in m^2."""
    h0_s_inv = h0_km_s_mpc_to_s_inv(h0_km_s_mpc)
    r_h = C_M_S / h0_s_inv
    return float(4.0 * PI * (r_h**2))


def hubble_entropy_log10(h0_km_s_mpc: float) -> float:
    """Return log10(S_H) where S_H = A_H/(4 ℓ_P^2)."""
    a_h = hubble_horizon_area_m2(h0_km_s_mpc)
    lp = planck_length_m()
    s_h = a_h / (4.0 * (lp**2))
    return math.log10(s_h)


@dataclass(frozen=True)
class SpinFoamVacuum:
    """Toy spin-foam-inspired vacuum energy mechanism.

    Parameters
    ----------
    gamma:
        Immirzi parameter (kept for future non-toy extensions).
    avg_amplitude:
        Dimensionless coarse-grained amplitude factor Ā.
    suppression_exponent_log10:
        Base-10 suppression exponent X in 10^{-X}. A value near ~122 is
        motivated by the horizon entropy scale.

    Notes
    -----
    This is a stand-in interface to make it easy to:
    - sweep X around 122–123,
    - quantify implied tuning in Ā, and
    - optionally use this as a target in constraint-optimization pipelines.
    """

    gamma: float = 0.2375
    avg_amplitude: float = 1.0
    suppression_exponent_log10: float = 122.0

    name: str = "spin_foam_vacuum"

    def __post_init__(self) -> None:
        if self.gamma <= 0:
            raise ValueError("gamma must be positive")
        if self.avg_amplitude <= 0:
            raise ValueError("avg_amplitude must be positive")
        if self.suppression_exponent_log10 < 0:
            raise ValueError("suppression_exponent_log10 must be >= 0")

    def describe_assumptions(self) -> str:
        return (
            "Toy spin-foam vacuum: ρ_Λ,eff = ρ_Pl * Ā * 10^{-X} with constant w=-1. "
            f"Parameters: gamma={self.gamma:.4g}, Ā={self.avg_amplitude:.4g}, X={self.suppression_exponent_log10:.4g}. "
            "This is not an EPRL/BC amplitude computation; it is a hierarchy/tuning probe."
        )

    def rho_lambda_eff_j_m3(self) -> float:
        rho_pl = planck_energy_density_j_m3()
        suppression = 10.0 ** (-self.suppression_exponent_log10)
        return float(rho_pl * self.avg_amplitude * suppression)

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:
        ensure_z_nonnegative(z)
        rho_de = self.rho_lambda_eff_j_m3()
        return MechanismOutput(
            result=MechanismResult(z=z, rho_de_j_m3=rho_de, w_de=-1.0),
            assumptions=self.describe_assumptions(),
        )
