from __future__ import annotations

import math
from dataclasses import dataclass

from ..constants import PI
from ..cosmology import h0_km_s_mpc_to_s_inv
from .base import CosmologyBackground, MechanismOutput, MechanismResult, ensure_z_nonnegative


@dataclass(frozen=True)
class RunningVacuumRVM:
    """Toy running-vacuum model (RVM) with a single parameter ν.

    This is scaffolding intended for constraint comparisons, not a derivation.

    A commonly used RVM ansatz is:

      ρ_Λ(H) = ρ_Λ0 + (3 ν / (8πG)) (H^2 - H0^2) * c^2

    Coupled self-consistently, one obtains a simple algebraic modification of H(z).
    Here we implement the standard toy closed form for flat backgrounds:

      H^2(z) = H0^2/(1-ν) * [Ωm (1+z)^3 + Ωr (1+z)^4 + ΩΛ - ν + Ωk (1+z)^2]

    and define ρ_DE(z) by assuming it equals Ω_DE(z) times the critical density.

    Notes:
    - For ν=0, this reduces to constant Λ (ρ_DE=ρ_Λ0).
    - This is a bookkeeping model; it does not enforce matter-DE exchange.
    """

    name: str = "running_vacuum_rvm"
    nu: float = 0.0

    def describe_assumptions(self) -> str:
        return (
            "Toy Running Vacuum Model: assumes an effective ν parameter that rescales H(z) via a simple algebraic form. "
            "Does not model explicit matter–DE exchange; intended for constraint scaffolding only."
        )

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:
        ensure_z_nonnegative(z)
        if not (-0.5 < self.nu < 0.5):
            raise ValueError("ν outside toy stability range")

        zp1 = 1.0 + z
        denom = 1.0 - self.nu
        if denom <= 0:
            raise ValueError("ν must be < 1 for this toy form")

        # E(z)^2 = (H/H0)^2
        ez2 = (
            bg.omega_m * (zp1**3)
            + bg.omega_r * (zp1**4)
            + bg.omega_k * (zp1**2)
            + bg.omega_lambda
            - self.nu
        ) / denom

        if ez2 <= 0:
            raise ValueError("Non-physical H(z)^2 encountered in toy RVM")

        # Define a toy evolving DE fraction by subtracting matter/radiation/curvature contributions.
        # This is not unique; we keep it simple and unit-consistent.
        omega_de_z = max(0.0, 1.0 - (bg.omega_m * (zp1**3) + bg.omega_r * (zp1**4) + bg.omega_k * (zp1**2)) / ez2)

        # ρ_c,E(z) scales as H(z)^2
        rho_de = bg.rho_lambda0_j_m3 * (omega_de_z / bg.omega_lambda) * ez2 if bg.omega_lambda > 0 else bg.rho_lambda0_j_m3

        # w is not well-defined without an explicit continuity equation; leave as None.
        return MechanismOutput(result=MechanismResult(z=z, rho_de_j_m3=rho_de, w_de=None), assumptions=self.describe_assumptions())
