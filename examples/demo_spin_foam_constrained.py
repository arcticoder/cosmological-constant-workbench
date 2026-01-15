#!/usr/bin/env python3
"""Demonstration: spin-foam-constrained mechanisms (predictor-independent path).

This demo shows how a toy spin-foam vacuum target can constrain free parameters
in phenomenological mechanisms without requiring the external LQG predictor.

Run:
  PYTHONPATH=src python examples/demo_spin_foam_constrained.py
"""

from ccw.integrations.lqg_constrained import (
    holographic_constrained_by_spin_foam,
    sequestering_constrained_by_spin_foam,
)
from ccw.mechanisms import CosmologyBackground
from ccw.mechanisms.spin_foam_vacuum import hubble_entropy_log10


def main() -> None:
    print("=" * 78)
    print("Spin-Foam-Constrained Mechanisms (predictor-independent)")
    print("=" * 78)
    print()

    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)

    log10_s = hubble_entropy_log10(bg.h0_km_s_mpc)
    print(f"H0 = {bg.h0_km_s_mpc:.2f} km/s/Mpc")
    print(f"Horizon entropy scale: log10(S_H) ≈ {log10_s:.3f}")
    print()

    print("Scenario 1: Holographic DE constrained by toy spin-foam vacuum")
    print("-" * 78)
    print()
    print("Target: ρ_Λ = ρ_Pl * Ā * 10^{-X} with X ≈ 122 (horizon entropy scale).")
    print()

    result_hde = holographic_constrained_by_spin_foam(
        bg,
        avg_amplitude=1.0,
        suppression_exponent_log10=122.0,
    )
    print(f"Spin-foam target: ρ_Λ = {result_hde.lqg_target_rho_j_m3:.3e} J/m³")
    print(f"Best c_factor:    {result_hde.best_fit_params['c_factor']:.2f}")
    print(f"Achieved:         ρ_DE = {result_hde.achieved_rho_j_m3:.3e} J/m³")
    print(f"Residual tuning:  log10(Δρ/ρ) = {result_hde.residual_tuning:.2f}")
    print(f"Success:          {result_hde.success}")
    print()
    print(result_hde.notes)
    print()
    print()

    print("Scenario 2: Sequestering constrained by toy spin-foam vacuum")
    print("-" * 78)
    print()
    print("Question: If spin-foam predicts ρ_Λ via horizon entropy, how much")
    print("cancellation f_cancel is required to match it from Planck-scale vacuum?")
    print()

    result_seq = sequestering_constrained_by_spin_foam(
        bg,
        avg_amplitude=1.0,
        suppression_exponent_log10=122.0,
        rho_vac_j_m3=1e113,
    )
    print(f"Spin-foam target: ρ_Λ = {result_seq.lqg_target_rho_j_m3:.3e} J/m³")
    print(f"Bare vacuum:      ρ_vac = 1.000e+113 J/m³")
    print(f"Required f_cancel: {result_seq.best_fit_params['f_cancel']:.15f}")
    print(f"Achieved:         ρ_DE = {result_seq.achieved_rho_j_m3:.3e} J/m³")
    print(f"Residual tuning:  ~10^{result_seq.residual_tuning:.0f}")
    print(f"Success:          {result_seq.success}")
    print()
    print(result_seq.notes)
    print()
    print()

    print("=" * 78)
    print("Conclusion:")
    print("=" * 78)
    print()
    print("- Holographic DE: Spin-foam target determines c_factor. If c ~ O(1),")
    print("  the combination is 'natural' modulo the spin-foam derivation itself.")
    print()
    print("- Sequestering: Extreme fine-tuning remains (~120 orders of magnitude).")
    print("  Spin-foam prediction does NOT eliminate the sequestering tuning problem.")
    print()
    print("Net assessment:")
    print("  This is a toy scaffolding (Phase K.27). The suppression exponent X and")
    print("  amplitude Ā are scan parameters, not derived from actual EPRL/BC vertex")
    print("  amplitudes. Next step: replace toy with real spin-foam amplitude calculation.")


if __name__ == "__main__":
    main()
