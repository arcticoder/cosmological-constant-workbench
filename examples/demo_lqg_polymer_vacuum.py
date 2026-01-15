#!/usr/bin/env python3
"""Demonstration: LQG polymer + toy derived vacuum term.

This extends the basic holonomy correction (which yields a negative effective
term -ρ^2/ρ_c) by adding a heuristic, constant-like positive contribution:

  ρ_Λ,LQG = ρ_Pl * prefactor * exp(-S_bounce/α)
  S_bounce = 4π * mu0_factor^2

This is NOT a first-principles derivation. It is a controlled way to:
- quantify how sensitive ρ_Λ,LQG is to the suppression parameters, and
- compute a tuning-like metric versus observed ρ_Λ,0.

Run:
  PYTHONPATH=src python examples/demo_lqg_polymer_vacuum.py
"""

from __future__ import annotations

import math

from ccw.mechanisms import CosmologyBackground
from ccw.mechanisms.lqg_polymer import (
    LQGPolymerDerivedVacuum,
    planck_energy_density_j_m3,
    toy_bounce_entropy,
    toy_derived_vacuum_energy_density_j_m3,
)


def log10_ratio(x: float, y: float) -> float:
    if y == 0:
        raise ZeroDivisionError("y must be nonzero")
    if x == 0:
        return float("-inf")
    return math.log10(abs(x / y))


def main() -> None:
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)

    print("=" * 78)
    print("LQG polymer + toy derived vacuum term")
    print("=" * 78)
    print(f"Observed ρΛ,0: {bg.rho_lambda0_j_m3:.3e} J/m³")
    print(f"Planck ρ_Pl:   {planck_energy_density_j_m3():.3e} J/m³")
    print()

    # A small grid of (mu0_factor, alpha) to show sensitivity.
    mu0_grid = [10.0, 30.0, 50.0, 60.0, 70.0]
    alpha_grid = [60.0, 120.0, 240.0]

    for alpha in alpha_grid:
        print("-" * 78)
        print(f"alpha_entropy = {alpha:.1f}")
        for mu0 in mu0_grid:
            s = toy_bounce_entropy(mu0_factor=mu0)
            rho_vac = toy_derived_vacuum_energy_density_j_m3(mu0_factor=mu0, alpha_entropy=alpha)

            mech = LQGPolymerDerivedVacuum(mu0_factor=mu0, alpha_entropy=alpha)
            rho_de0 = mech.effective_rho_de_j_m3(0.0, bg)

            metric = log10_ratio(rho_vac, bg.rho_lambda0_j_m3)
            print(
                f"  mu0={mu0:>5.1f}  S={s:>9.2e}  ρ_vac={rho_vac:>10.3e}  "
                f"log10(ρ_vac/ρΛ,0)={metric:>7.2f}  ρ_DE,total(z=0)={rho_de0:>10.3e}"
            )

    print("-" * 78)
    print("Interpretation:")
    print("- The holonomy correction alone does not generate positive Λ.")
    print("- The derived-vacuum toy can hit ρΛ,0 only for specific (mu0, α) choices,")
    print("  so it should be treated as a tuning quantifier, not a prediction.")


if __name__ == "__main__":
    main()
