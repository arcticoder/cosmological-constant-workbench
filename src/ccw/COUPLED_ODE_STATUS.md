# Coupled ODE Solver Status

## Current State (Jan 13, 2026)

The coupled ODE solver (`src/ccw/coupled_ode.py`) is **NOT WORKING** due to fundamental unit/normalization issues.

### Problem

The solver attempts to integrate field equations (φ, Π) coupled to the Friedmann equation for H(a), but encounters:

1. **Numerical stiffness**: Integration requires step sizes smaller than machine precision
2. **Unit mismatch**: Mixing dimensionless field variables (φ in Planck units) with SI energy densities
3. **Overflow issues**: Exponential potential exp(-λφ) overflows when φ → -∞ during backward integration

### Root Cause

The field φ and its conjugate momentum Π are dimensionless (normalized to M_Pl), but the potential V(φ) and derivative V'(φ) must be converted to SI energy densities (J/m³) for use in the Friedmann equation. The current implementation:

```python
rho_phi = 0.5 * Pi**2 + V_phi  # Dimensionless
rho_tot = rho_m + rho_r + rho_phi * (C_M_S**2)  # Mixing units!
```

This creates a huge scale mismatch: `C_M_S**2 ~ 9×10^16`, causing the field energy to dominate immediately and the integrator to fail.

### What Needs to be Fixed

1. **Proper normalization**: Define φ and Π in units such that V(φ) directly gives ρ_DE in J/m³
   - Option A: φ in Planck masses, V(φ) = V₀ M_Pl⁴ → scale by M_Pl⁴c² → J/m³
   - Option B: φ dimensionless, V(φ) = V₀ ρ_crit → already in J/m³

2. **Consistent derivatives**: V'(φ) units must match (∂V/∂φ) × ρ_crit for field equation

3. **Initial conditions**: φ₀, Π₀ must be small enough to avoid runaway but large enough to matter

4. **Integration method**: May need adaptive stiffness handling or event detection for φ → 0

### Recommended Path Forward

**Short term** (to unblock progress):
- Skip the coupled ODE solver for now
- Continue using `MechanismHz` with algebraic H(z) from ρ_DE(z) on fixed background
- Note in `ccw-report` that self-consistency is approximate

**Medium term** (after emergent/GW/backreaction are working):
- Revisit coupled ODE with proper dimensional analysis
- Benchmark against known quintessence solutions (e.g., Wetterich exponential)
- Use dimensionless variables throughout: τ = ln(a), x = φ/M_Pl, etc.

**Long term** (if needed for discovery):
- Implement event-driven integration to detect instabilities
- Add automated parameter search to find viable (φ₀, Π₀) pairs
- Compare to CLASS/CAMB quintessence modules for validation

## Files Affected

- [x] `src/ccw/coupled_ode.py` (created but not functional)
- [x] `tests/test_coupled_ode.py` (8/10 tests failing)
- [ ] `src/ccw/mechanisms/scalar_field.py` (integration pending)
- [ ] `examples/demo_coupled_ode.py` (not created)

## Next Steps

1. **Mark J.22 as "blocked" in TASKS.md**
2. **Move to H.18 (emergent gravity)** which doesn't require coupled ODE
3. **Re-evaluate J.22 after H.18, I.21, J.23 are done**

---

**Commit message**: "WIP: Coupled ODE solver - unit normalization issues, tests failing (8/10)"

**Lessons learned**: Always verify dimensional analysis before implementing physics equations. Planck unit conventions require explicit conversion factors when mixing with SI.
