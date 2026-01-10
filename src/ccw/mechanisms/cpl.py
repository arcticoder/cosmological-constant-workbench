from __future__ import annotations

import math
from dataclasses import dataclass

from .base import CosmologyBackground, MechanismOutput, MechanismResult, ensure_z_nonnegative


@dataclass(frozen=True)
class CPLQuintessence:
    """CPL parameterization as a stand-in for quintessence-like dynamics.

    This is *not* a derived scalar-field solution; it is a commonly used
    phenomenological form:

      w(a) = w0 + wa (1 - a), where a = 1/(1+z)

    The implied density evolution is:

      ρ_DE(a) = ρ0 * a^{-3(1+w0+wa)} * exp(3 wa (a - 1))

    References: Chevallier–Polarski–Linder (CPL).
    """

    name: str = "cpl_quintessence"
    w0: float = -1.0
    wa: float = 0.0

    def describe_assumptions(self) -> str:
        return (
            "Phenomenological CPL dark-energy parameterization (not a specific scalar-field potential). "
            "Assumes smooth DE with equation-of-state w(a)=w0+wa(1-a)."
        )

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:
        ensure_z_nonnegative(z)
        a = 1.0 / (1.0 + z)

        rho0 = bg.rho_lambda0_j_m3
        exponent = -3.0 * (1.0 + self.w0 + self.wa)
        rho = rho0 * (a**exponent) * math.exp(3.0 * self.wa * (a - 1.0))

        w = self.w0 + self.wa * (1.0 - a)
        return MechanismOutput(result=MechanismResult(z=z, rho_de_j_m3=rho, w_de=w), assumptions=self.describe_assumptions())
