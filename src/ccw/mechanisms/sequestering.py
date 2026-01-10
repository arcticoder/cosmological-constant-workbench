from __future__ import annotations

from dataclasses import dataclass

from .base import CosmologyBackground, MechanismOutput, MechanismResult, ensure_z_nonnegative


@dataclass(frozen=True)
class SequesteringToy:
    """Vacuum energy sequestering toy model.

    Real sequestering proposals modify how vacuum energy gravitates (global
    constraints / auxiliary fields). Here we only model an *effective* residual
    after cancellation:

      ρ_DE(z) = ρ_Λ0 + δρ

    with δρ as a user-controlled parameter.
    """

    name: str = "sequestering_toy"
    delta_rho_j_m3: float = 0.0

    def describe_assumptions(self) -> str:
        return (
            "Toy sequestering: assumes vacuum contributions cancel leaving a small residual δρ added to ρ_Λ0. "
            "No underlying dynamics implemented."
        )

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:
        ensure_z_nonnegative(z)
        rho = bg.rho_lambda0_j_m3 + self.delta_rho_j_m3
        return MechanismOutput(result=MechanismResult(z=z, rho_de_j_m3=rho, w_de=-1.0), assumptions=self.describe_assumptions())
