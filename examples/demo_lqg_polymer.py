#!/usr/bin/env python3
"""Demonstration: LQG polymer (holonomy) corrections as effective dark energy.

This example is intentionally conservative:
- It does NOT claim a cosmological-constant solution.
- It shows how large a rescaling of the polymer scale (μ0) would be required to
  make holonomy corrections visible at late times.

We treat the LQC correction
  ρ_DE,eff(z) = -ρ(z)^2 / ρ_c
as an *effective* dark-energy density within the workbench mechanism interface.

Run:
  PYTHONPATH=src python examples/demo_lqg_polymer.py
"""

from __future__ import annotations

import math

from ccw.mechanisms import CosmologyBackground
from ccw.mechanisms.lqg_polymer import LQGPolymerCosmology, planck_energy_density_j_m3


def main() -> None:
    bg = CosmologyBackground(h0_km_s_mpc=67.4, omega_lambda=0.6889, omega_m=0.3111, omega_r=0.0)

    print("=" * 78)
    print("LQG polymer cosmology: effective ρ_DE,eff(z) = -ρ(z)^2 / ρ_c")
    print("=" * 78)
    print(f"Observed ρΛ,0: {bg.rho_lambda0_j_m3:.3e} J/m³")
    print(f"Planck ρ_Pl:   {planck_energy_density_j_m3():.3e} J/m³")
    print()

    z_values = [0.0, 1.0, 10.0, 1000.0]
    mu0_factors = [1.0, 1e30, 1e50, 1e60, 1e61]

    for mu0 in mu0_factors:
        mech = LQGPolymerCosmology(mu0_factor=float(mu0))
        rho_c = mech.rho_c_j_m3()

        print("-" * 78)
        print(f"μ0_factor = {mu0:.1e}  ->  ρ_c = {rho_c:.3e} J/m³")
        for z in z_values:
            rho_de = mech.effective_rho_de_j_m3(float(z), bg)
            ratio = rho_de / bg.rho_lambda0_j_m3
            sign = "+" if rho_de >= 0 else "-"
            print(
                f"  z={z:<6g}  ρ_DE,eff={rho_de:.3e} J/m³  "
                f"(ρ_DE,eff/ρΛ,0={ratio:.3e}, sign={sign})"
            )

    print("-" * 78)
    print("Notes:")
    print("- For μ0_factor≈1 (Planck-scale), ρ_DE,eff at late times is negligible.")
    print("- Making μ0_factor enormous lowers ρ_c and amplifies corrections, but")
    print("  the effective density here is negative and strongly redshift-dependent.")
    print("- A positive late-time Λ would require additional physics/terms beyond this")
    print("  minimal holonomy correction bookkeeping.")


if __name__ == "__main__":
    main()
