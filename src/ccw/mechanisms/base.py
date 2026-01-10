from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Protocol

from ..cosmology import observed_lambda_from_h0_omega


def ensure_z_nonnegative(z: float) -> None:
    if z < 0:
        raise ValueError("Redshift z must be >= 0")


@dataclass(frozen=True)
class CosmologyBackground:
    """Minimal background needed for toy dark-energy mechanisms.

    This is deliberately lightweight: it supports mapping from (H0, ΩΛ) to
    observed ρ_Λ today, and provides Ωm/Ωr/Ωk for simple H(z) bookkeeping.
    """

    h0_km_s_mpc: float = 67.4
    omega_lambda: float = 0.6889
    omega_m: float = 0.3111
    omega_r: float = 0.0
    omega_k: float = 0.0

    @property
    def rho_lambda0_j_m3(self) -> float:
        return observed_lambda_from_h0_omega(self.h0_km_s_mpc, self.omega_lambda).rho_lambda_j_m3


@dataclass(frozen=True)
class MechanismResult:
    z: float
    rho_de_j_m3: float
    w_de: Optional[float]


@dataclass(frozen=True)
class MechanismOutput:
    """For future extension: arrays, extra observables, diagnostics."""

    result: MechanismResult
    assumptions: str


class Mechanism(Protocol):
    name: str

    def describe_assumptions(self) -> str:  # pragma: no cover
        ...

    def evaluate(self, z: float, bg: CosmologyBackground) -> MechanismOutput:  # pragma: no cover
        ...
