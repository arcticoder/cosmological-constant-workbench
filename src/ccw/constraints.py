from __future__ import annotations

from dataclasses import dataclass

from .constants import C_M_S, G_M3_KG_S2, PI
from .cosmology import ObservedLambdaResult, observed_lambda_from_h0_omega


@dataclass(frozen=True)
class LambdaObservables:
    """Convenience container for (ΩΛ, H0) <-> (ρΛ, Λ) conversions."""

    h0_km_s_mpc: float
    omega_lambda: float
    rho_lambda_j_m3: float
    lambda_m_inv2: float


def from_h0_omega(h0_km_s_mpc: float, omega_lambda: float) -> LambdaObservables:
    obs: ObservedLambdaResult = observed_lambda_from_h0_omega(h0_km_s_mpc, omega_lambda)
    return LambdaObservables(
        h0_km_s_mpc=h0_km_s_mpc,
        omega_lambda=omega_lambda,
        rho_lambda_j_m3=obs.rho_lambda_j_m3,
        lambda_m_inv2=obs.lambda_m_inv2,
    )


def omega_lambda_from_lambda_h0(lambda_m_inv2: float, h0_s_inv: float) -> float:
    """ΩΛ = Λ c^2 / (3 H0^2)."""

    if lambda_m_inv2 < 0:
        raise ValueError("Λ must be >= 0 for this bookkeeping conversion")
    if h0_s_inv <= 0:
        raise ValueError("H0 must be positive")
    return lambda_m_inv2 * (C_M_S**2) / (3.0 * (h0_s_inv**2))


def lambda_from_rho_lambda(rho_lambda_j_m3: float) -> float:
    """Λ = 8πG ρ / c^4 when ρ is an energy density (J/m^3)."""

    if rho_lambda_j_m3 < 0:
        raise ValueError("ρ_Λ must be >= 0 for this bookkeeping conversion")
    return (8.0 * PI * G_M3_KG_S2 * rho_lambda_j_m3) / (C_M_S**4)


def rho_lambda_from_lambda(lambda_m_inv2: float) -> float:
    """ρ = Λ c^4 / (8πG) when ρ is an energy density (J/m^3)."""

    if lambda_m_inv2 < 0:
        raise ValueError("Λ must be >= 0 for this bookkeeping conversion")
    return (lambda_m_inv2 * (C_M_S**4)) / (8.0 * PI * G_M3_KG_S2)


def sanity_check_background(omega_m: float, omega_lambda: float, omega_r: float = 0.0, omega_k: float = 0.0) -> None:
    """Basic bounds and consistency checks (not a cosmological fit)."""

    for name, val in [("Ωm", omega_m), ("ΩΛ", omega_lambda), ("Ωr", omega_r), ("Ωk", omega_k)]:
        if not (-5.0 <= val <= 5.0):
            raise ValueError(f"{name} outside sanity range")

    total = omega_m + omega_lambda + omega_r + omega_k
    if not (0.0 <= total <= 5.0):
        raise ValueError("Ω total outside sanity range")
