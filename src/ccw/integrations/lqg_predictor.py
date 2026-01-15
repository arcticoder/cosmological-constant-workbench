"""Optional adapter for comparing against lqg-cosmological-constant-predictor.

This adapter is guarded: all imports are lazy and all functions handle import
failures gracefully, returning availability metadata.

Usage:
    from ccw.integrations.lqg_predictor import lqg_predictor_available, compare_baseline

    if lqg_predictor_available():
        result = compare_baseline(h0, omega_lambda)
        print(result)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .base import try_import_lqg_cosmological_constant_predictor


@dataclass(frozen=True)
class LQGPredictorComparison:
    """Comparison result between CCW baseline and LQG predictor."""

    lqg_lambda_m_minus2: float
    lqg_rho_j_m3: float
    ccw_lambda_m_minus2: float
    ccw_rho_j_m3: float
    ratio_lambda: float
    ratio_rho: float
    notes: str


@dataclass(frozen=True)
class LQGPredictorFirstPrinciplesConfig:
    """Configuration for running the external LQG predictor.

    Notes
    -----
    This is intentionally a thin wrapper around the external API. We expose a
    generic `params_overrides` dict so we can sweep external parameters without
    pinning CCW to the predictor's internal dataclass layout.
    """

    target_scale_m: float = 1e-15
    include_uncertainty: bool = False
    params_overrides: Optional[dict[str, Any]] = None


@dataclass(frozen=True)
class LQGPredictorFirstPrinciplesResult:
    """Stable CCW-side summary of an external LQG predictor run."""

    lambda_effective_m_minus2: float
    vacuum_energy_density_j_m3: float
    mu_scale: Optional[float]
    enhancement_factor: Optional[float]
    scale_correction: Optional[float]
    notes: str


def _apply_params_overrides(params: Any, overrides: dict[str, Any]) -> None:
    for key, value in overrides.items():
        if not hasattr(params, key):
            raise ValueError(f"Unknown predictor param override: {key}")
        setattr(params, key, value)


def run_first_principles(
    config: LQGPredictorFirstPrinciplesConfig = LQGPredictorFirstPrinciplesConfig(),
) -> LQGPredictorFirstPrinciplesResult:
    """Run the external predictor and return a stable CCW-side summary.

    Raises
    ------
    ImportError
        If the external repo is not importable.
    ValueError
        If params_overrides contain unknown keys.
    """

    if not lqg_predictor_available():
        raise ImportError(
            "lqg-cosmological-constant-predictor not available. "
            "Ensure it is installed or in the Python path."
        )

    from cosmological_constant_predictor import CosmologicalConstantPredictor, CosmologicalParameters

    params = CosmologicalParameters()
    if config.params_overrides:
        _apply_params_overrides(params, config.params_overrides)

    predictor = CosmologicalConstantPredictor(params)
    result = predictor.predict_lambda_from_first_principles(
        target_scale=config.target_scale_m,
        include_uncertainty=config.include_uncertainty,
    )

    # Extract fields defensively; external API may add/remove attributes.
    mu_scale = getattr(result, "mu_scale", None)
    enhancement_factor = getattr(result, "enhancement_factor", None)
    scale_correction = getattr(result, "scale_correction", None)

    notes = (
        "External LQG predictor run: vacuum_energy_density comes from the predictor's "
        "enhanced polymer-modified vacuum-energy sum with SU(2) corrections and scale-dependent Λ_eff(ℓ)."
    )

    return LQGPredictorFirstPrinciplesResult(
        lambda_effective_m_minus2=float(result.lambda_effective),
        vacuum_energy_density_j_m3=float(result.vacuum_energy_density),
        mu_scale=None if mu_scale is None else float(mu_scale),
        enhancement_factor=None if enhancement_factor is None else float(enhancement_factor),
        scale_correction=None if scale_correction is None else float(scale_correction),
        notes=notes,
    )


def lqg_predictor_available() -> bool:
    """Check if lqg-cosmological-constant-predictor is importable."""

    return try_import_lqg_cosmological_constant_predictor().available


def compare_baseline(
    h0_km_s_mpc: float = 67.4, omega_lambda: float = 0.6889
) -> LQGPredictorComparison:
    """Compare CCW baseline (H0, ΩΛ) conversion with LQG predictor first-principles result.

    Args:
        h0_km_s_mpc: Hubble constant in km/s/Mpc
        omega_lambda: Dark-energy density fraction today

    Returns:
        Comparison result.

    Raises:
        ImportError: if lqg-cosmological-constant-predictor is not available.
    """

    if not lqg_predictor_available():
        raise ImportError(
            "lqg-cosmological-constant-predictor not available. "
            "Ensure it is installed or in the Python path."
        )

    # Import CCW baseline
    from ..cosmology import observed_lambda_from_h0_omega

    # Lazy import of LQG predictor
    from cosmological_constant_predictor import CosmologicalConstantPredictor

    # CCW baseline (observational conversion)
    ccw_result = observed_lambda_from_h0_omega(h0_km_s_mpc, omega_lambda)

    # LQG predictor first-principles result
    predictor = CosmologicalConstantPredictor()
    # Note: The LQG predictor returns a complex result object; we extract the key fields.
    # This may require adjustment if the LQG API changes.
    lqg_full_result = predictor.predict_lambda_from_first_principles(
        include_uncertainty=False
    )

    # Extract LQG predictions
    lqg_lambda = lqg_full_result.lambda_effective
    lqg_rho = lqg_full_result.vacuum_energy_density

    # Compute ratios
    ratio_lambda = lqg_lambda / ccw_result.lambda_m_inv2
    ratio_rho = lqg_rho / ccw_result.rho_lambda_j_m3

    notes = (
        "LQG predictor uses first-principles polymer quantization + SU(2) 3nj corrections; "
        "CCW baseline uses standard ΛCDM conversion from (H0, ΩΛ)."
    )

    return LQGPredictorComparison(
        lqg_lambda_m_minus2=lqg_lambda,
        lqg_rho_j_m3=lqg_rho,
        ccw_lambda_m_minus2=ccw_result.lambda_m_inv2,
        ccw_rho_j_m3=ccw_result.rho_lambda_j_m3,
        ratio_lambda=ratio_lambda,
        ratio_rho=ratio_rho,
        notes=notes,
    )
