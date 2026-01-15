#!/usr/bin/env python3
"""Demonstration: toy spin-foam vacuum mechanism (Phase K.27 scaffolding).

Run:
  PYTHONPATH=src python examples/demo_spin_foam_vacuum.py

This is a toy parameterization, not an EPRL/BC amplitude computation.
"""

import math

from ccw.mechanisms import CosmologyBackground, SpinFoamVacuum
from ccw.mechanisms.lqg_polymer import planck_energy_density_j_m3
from ccw.mechanisms.spin_foam_vacuum import hubble_entropy_log10


def main() -> None:
    print("=" * 78)
    print("Toy Spin-Foam Vacuum: entropy-motivated suppression scale")
    print("=" * 78)
    print()

    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111)

    rho_obs = bg.rho_lambda0_j_m3
    log10_s = hubble_entropy_log10(bg.h0_km_s_mpc)

    rho_pl = planck_energy_density_j_m3()
    print(f"H0 = {bg.h0_km_s_mpc:.2f} km/s/Mpc")
    print(f"log10(S_H) ≈ {log10_s:.3f}")
    print(f"log10(ρ_Pl [J/m³]) = {math.log10(rho_pl):.3f}")
    print(f"Observed ρ_Λ,0 = {rho_obs:.3e} J/m³")
    print()

    print("Scan suppression exponent X in ρ = ρ_Pl * Ā * 10^{-X} (Ā=1)")
    print("-" * 78)
    for x in [120.0, 121.0, 122.0, 122.5, 123.0, 124.0]:
        mech = SpinFoamVacuum(avg_amplitude=1.0, suppression_exponent_log10=x)
        rho = mech.rho_lambda_eff_j_m3()
        ratio = rho / rho_obs if rho_obs > 0 else float("nan")
        print(f"X={x:5.1f}  ρ={rho: .3e} J/m³   ρ/ρ_obs={ratio: .3e}")

    print()
    print("Interpretation:")
    print("- If X is 'derived' (e.g., tied to horizon entropy), then Ā~O(1) matching")
    print("  the observed scale would look less like independent fine-tuning.")
    print("- In this toy, both X and Ā are scan knobs; K.27’s next step is replacing")
    print("  them with an actual amplitude calculation or a defensible proxy.")


if __name__ == "__main__":
    main()
