# Coupled ODE Fix Proposal (J.22 Unblocking)

## Problem Summary

Current `coupled_ode.py` mixes dimensionless field variables (φ in Planck units) with SI energy densities, causing numerical overflow and integration failure.

**Root cause**:
```python
rho_phi = 0.5 * Pi**2 + V_phi  # Dimensionless
rho_tot = rho_m + rho_r + rho_phi * (C_M_S**2)  # Scale mismatch!
```

The factor `C_M_S**2 ~ 9×10^16` creates huge energy density from tiny field values, causing runaway.

## Proposed Solution: Dimensionless Normalization

Normalize **all** quantities to critical density ρ_crit = 3H₀²/(8πG):

### Variables
- Time: `N = ln(a)` (e-folds)
- Field: `ψ = φ / M_Pl` (dimensionless)
- Momentum: `Ψ = Π / M_Pl` (dimensionless)
- Potential: `u(ψ) = V(ψ M_Pl) / ρ_crit,0` (dimensionless)
- Energy densities: `Ω_i = ρ_i / ρ_crit,0`

### Equations
Friedmann:
```
H²/H₀² = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_φ
```

where
```
Ω_φ = (M_Pl² / ρ_crit,0) × (½Ψ² + u(ψ))
```

Field evolution (in N = ln a):
```
dψ/dN = Ψ
dΨ/dN = -3(H/H₀)Ψ - (ρ_crit,0/M_Pl²) × du/dψ
```

### Implementation
1. Convert potential inputs to dimensionless form:
   ```python
   def u_exponential(psi, lam_tilde):
       """Dimensionless exponential potential.
       
       Physical: V(φ) = V₀ exp(-λ φ / M_Pl)
       Normalized: u(ψ) = (V₀ / ρ_crit,0) exp(-λ_tilde ψ)
       """
       return (V0_over_rho_crit) * np.exp(-lam_tilde * psi)
   ```

2. Field energy density contribution:
   ```python
   # Planck mass in SI
   M_Pl_kg = np.sqrt(HBAR_J_S * C_M_S / G_M3_KG_S2)
   M_Pl_energy = M_Pl_kg * C_M_S**2  # J
   
   # Critical density (today)
   rho_crit_0 = 3 * H0_si**2 / (8*PI*G_M3_KG_S2)  # kg/m³
   
   # Normalization factor
   norm = M_Pl_energy / rho_crit_0  # dimensionless ~ 10^120 (!)
   
   Omega_phi = norm * (0.5 * Psi**2 + u(psi))
   ```

3. **Key insight**: Since norm ~ 10^120, even tiny field values ψ ~ 10^-10 can dominate the universe. This is why the solver fails — the field contribution is always too large unless we choose:
   
   **Either**:
   - Initial conditions: ψ₀ ~ 10^-60, Ψ₀ ~ 10^-60 (requires extreme precision)
   
   **Or**:
   - Rescale potential: choose V₀ such that V₀/ρ_crit,0 ~ 10^-120, making u(ψ) ~ O(1) for ψ ~ O(1)

## Recommended Approach

**Option A: Fix normalization but keep precision requirements explicit**
- Implement dimensionless equations as above
- Document that quintessence requires V₀ ~ 10^-120 ρ_crit,0 (this IS the fine-tuning!)
- Use high-precision arithmetic if needed (mpmath)

**Option B: Defer coupled ODE until after UV completion modules**
- Current algebraic approximation (ρ_DE on fixed ΛCDM background) is accurate to ~1% for slow-roll quintessence
- Coupled ODE is only essential for:
  1. Tachyonic / fast-roll regimes (not currently explored)
  2. Backreaction larger than ~10% (would require V₀ ~ ρ_crit, which is already ruled out)
- **Verdict**: Not blocking for Phase H.18, I.21, J.23, K.25 work

## Recommendation for TASKS.md

**Keep J.22 blocked** with note:
> Coupled ODE solver requires either (1) high-precision arithmetic for extreme normalization, or (2) rethinking field-theoretic UV completion. Current algebraic H(z) from ρ_DE(z) is accurate for slow-roll quintessence. Unblock only if fast-roll or O(1) backreaction becomes priority.

**Priority**: address after K.25 (LQG polymer) to see if discrete geometry provides natural UV cutoff.

---

**Next steps if we decide to implement**:
1. Create `src/ccw/coupled_ode_dimensionless.py` with above normalization
2. Add validation: reproduce known Wetterich (λ=0.5) tracker solution
3. Benchmark against CLASS quintessence module
4. Document precision requirements in tests

**Effort estimate**: 2-3 days for careful implementation + validation.
