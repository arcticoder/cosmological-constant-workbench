from __future__ import annotations

from dataclasses import dataclass

from .base import CosmologyBackground, MechanismOutput, MechanismResult, ensure_z_nonnegative


@dataclass(frozen=True)
class UnimodularBookkeeping:
    """Unimodular gravity bookkeeping stub.

    In unimodular gravity, Λ can appear as an integration constant rather than a
    parameter in the action. This stub does *not* implement a quantum mechanism
    fixing that constant; it only provides a placeholder handle for later work.

    Current behavior: returns constant ρ_DE = ρ_Λ0.
    """

    name: str = "unimodular_bookkeeping"

    def describe_assumptions(self) -> str:
        return (
            "Bookkeeping stub for unimodular gravity: treats Λ as an integration constant. "
            "No dynamical selection mechanism implemented; returns constant ρ_DE today."
        )

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:
        ensure_z_nonnegative(z)
        return MechanismOutput(
            result=MechanismResult(z=z, rho_de_j_m3=bg.rho_lambda0_j_m3, w_de=-1.0),
            assumptions=self.describe_assumptions(),
        )
