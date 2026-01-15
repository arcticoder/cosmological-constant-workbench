#!/usr/bin/env python3
"""Demonstration: constrained mechanisms using lqg-cosmological-constant-predictor.

This demo uses the external predictor's first-principles vacuum-energy output
as the target ρ_Λ to constrain phenomenological mechanisms.

Run:
  PYTHONPATH=src python examples/demo_lqg_predictor_constrained.py
"""

from ccw.integrations.lqg_constrained import (
    holographic_constrained_by_lqg_predictor,
    sequestering_constrained_by_lqg_predictor,
)
from ccw.integrations.lqg_predictor import lqg_predictor_available, run_first_principles
from ccw.integrations.lqg_predictor import LQGPredictorFirstPrinciplesConfig
from ccw.mechanisms import CosmologyBackground


def main() -> None:
    print("=" * 78)
    print("LQG Predictor-Constrained Mechanisms")
    print("=" * 78)
    print()

    if not lqg_predictor_available():
        print("⚠ LQG predictor not available.")
        print("This demo requires lqg-cosmological-constant-predictor to be importable.")
        return

    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)

    cfg = LQGPredictorFirstPrinciplesConfig(
        target_scale_m=1e-15,
        include_uncertainty=False,
        # Example override knobs (optional):
        # params_overrides={"mu_polymer": 0.15, "gamma_coefficient": 1.0, "volume_eigenvalue_cutoff": 10.0},
        params_overrides=None,
    )

    fp = run_first_principles(cfg)
    print(f"Predictor target scale: {cfg.target_scale_m:.2e} m")
    print(f"Predictor Λ_eff:        {fp.lambda_effective_m_minus2:.3e} m⁻²")
    print(f"Predictor ρ_vac:        {fp.vacuum_energy_density_j_m3:.3e} J/m³")
    if fp.mu_scale is not None:
        print(f"Predictor μ(scale):     {fp.mu_scale:.4g}")
    print()

    print("Scenario 1: Holographic DE constrained by predictor")
    print("-" * 78)
    res_h = holographic_constrained_by_lqg_predictor(bg, target_scale_m=cfg.target_scale_m)
    print(f"Target:   ρ_Λ = {res_h.lqg_target_rho_j_m3:.3e} J/m³")
    print(f"Best c:   {res_h.best_fit_params.get('c_factor', float('nan')):.3g}")
    print(f"Achieved: ρ_DE = {res_h.achieved_rho_j_m3:.3e} J/m³")
    print(f"Tuning:   log10(Δρ/ρ) = {res_h.residual_tuning:.2f}")
    print(f"Success:  {res_h.success}")
    print(res_h.notes)
    print()

    print("Scenario 2: Sequestering constrained by predictor")
    print("-" * 78)
    res_s = sequestering_constrained_by_lqg_predictor(bg, target_scale_m=cfg.target_scale_m, rho_vac_j_m3=1e113)
    print(f"Target:   ρ_Λ = {res_s.lqg_target_rho_j_m3:.3e} J/m³")
    print(f"f_cancel: {res_s.best_fit_params.get('f_cancel', float('nan')):.15f}")
    print(f"Tuning:   ~10^{res_s.residual_tuning:.0f}")
    print(f"Success:  {res_s.success}")
    print(res_s.notes)


if __name__ == "__main__":
    main()
