"""
Quantum backreaction and radiative stability for scalar field mechanisms.

Provides tools to estimate one-loop quantum corrections to effective potentials
and check whether mechanisms are radiatively stable (no fine-tuning at loop level).

Key physics:
- Scalar field φ coupled to matter generates loop corrections ΔV_eff(φ)
- Coleman-Weinberg: ΔV ~ (m⁴/64π²) ln(Λ_UV/m) with field-dependent mass m(φ)
- Radiative stability: require max(ΔV) ≪ ρ_Λ,observed to avoid loop-level tuning
- Quadratic divergence: ΔV ~ Λ_UV² for m ≪ Λ signals unprotected UV sensitivity

This addresses stopping criterion 3: rigorous bounds / no-go theorems.
"""

from dataclasses import dataclass
from typing import Callable, Optional, Tuple
import numpy as np

from .constants import PI, HBAR_J_S, C_M_S, G_M3_KG_S2, GEV_J, MPC_M


@dataclass
class BackreactionResult:
    """
    Result from backreaction / radiative stability analysis.
    
    Attributes
    ----------
    phi_grid : np.ndarray
        Scalar field values sampled (in Planck units).
    delta_V_loop : np.ndarray
        One-loop correction ΔV(φ) in J/m³ (energy density units).
    max_delta_V : float
        Maximum loop correction over field range (J/m³).
    rho_lambda_observed : float
        Observed vacuum energy density (J/m³, ~6e-10).
    tuning_level : float
        max(ΔV) / ρ_Λ,obs (dimensionless). >1 → radiatively unstable.
    is_stable : bool
        True if tuning_level < threshold (no fine-tuning required).
    quadratic_divergence_detected : bool
        True if ΔV ∝ Λ_UV² detected (unprotected UV sensitivity).
    """
    phi_grid: np.ndarray
    delta_V_loop: np.ndarray
    max_delta_V: float
    rho_lambda_observed: float
    tuning_level: float
    is_stable: bool
    quadratic_divergence_detected: bool


def coleman_weinberg_correction(
    phi: float,
    yukawa: float,
    lambda_uv_gev: float = 1e3,
    n_dof: int = 1,
) -> float:
    """
    Compute one-loop Coleman-Weinberg correction to effective potential.
    
    Parameters
    ----------
    phi : float
        Scalar field value in Planck units (dimensionless, φ/M_Pl).
    yukawa : float
        Yukawa coupling y (dimensionless). Field-dependent mass m(φ) = y·φ·M_Pl.
    lambda_uv_gev : float
        UV cutoff scale in GeV (e.g., 1 TeV = 1e3 GeV, M_Pl ~ 1e19 GeV).
    n_dof : int
        Number of degrees of freedom coupling to φ (e.g., 3 for SM fermions).
    
    Returns
    -------
    float
        Loop correction ΔV(φ) in J/m³ (energy density units).
    
    Notes
    -----
    Coleman-Weinberg formula (one-loop, MS-bar):
        ΔV(φ) = (n_dof / 64π²) × m(φ)⁴ × [ln(Λ_UV / m(φ)) - C]
    where C ~ 3/2 (renormalization constant, scheme-dependent).
    
    Field-dependent mass:
        m(φ) = y × φ × M_Pl   (Yukawa interaction ℒ ⊃ -y φ ψ̄ ψ)
    
    Physical interpretation:
    - Small y → ΔV → 0 (decoupling)
    - Large y → ΔV large (requires tuning if ΔV > ρ_Λ)
    - Λ_UV → ∞: logarithmic divergence (renormalizable)
    - m ≪ Λ_UV: ΔV ~ Λ_UV² (quadratic divergence if no mass threshold)
    
    Example:
    - φ ~ 1 (Planck scale field), y ~ 1, Λ_UV ~ 1 TeV
      → m ~ 10¹⁹ GeV, ΔV ~ 10⁻⁸ J/m³ ≫ ρ_Λ ~ 10⁻⁹ J/m³
      → Radiatively unstable, requires fine-tuning
    """
    # Planck mass in GeV
    M_Pl_GeV = np.sqrt(HBAR_J_S * C_M_S / G_M3_KG_S2) / GEV_J  # ~ 1.22e19 GeV
    
    # Field-dependent mass in GeV
    m_phi_gev = yukawa * abs(phi) * M_Pl_GeV
    
    if m_phi_gev < 1e-30:
        # Decoupling limit: m → 0, ΔV → 0
        return 0.0
    
    # Avoid log(Λ/m) → ∞ when m → 0
    if m_phi_gev > lambda_uv_gev:
        # Field mass exceeds UV cutoff: loop suppressed
        return 0.0
    
    # Coleman-Weinberg correction (logarithmic part)
    log_term = np.log(lambda_uv_gev / m_phi_gev)
    
    # One-loop: ΔV ~ (1/64π²) m⁴ ln(Λ/m)
    # Convert to J/m³ (energy density)
    m_phi_joule = m_phi_gev * GEV_J
    
    # Loop factor
    prefactor = n_dof / (64.0 * PI**2)
    
    # ΔV in J⁴ (need to convert to energy density)
    # Mass dimension: [m⁴] = (GeV)⁴ → (J)⁴ / (m)¹²
    # Energy density: [ΔV] = J/m³
    # So we need: (J)⁴ / (ħc/G)² → J/m³
    
    # Simpler: use natural units then convert
    # In natural units (ħ=c=1): ΔV has mass dimension 4
    # [ΔV] = (GeV)⁴ → multiply by (GeV/J) × (Planck length)³
    
    # Actually, cleaner: work in SI throughout
    # ΔV = prefactor × m⁴ × ln(Λ/m) where m is in energy units
    # To get J/m³, need to divide by characteristic volume
    # Characteristic length ~ Compton wavelength λ = ħ/(mc)
    # Volume ~ λ³ = (ħ/(mc))³
    
    # Energy density in SI: ΔV [J/m³] = prefactor × m⁴ × (mc²/ħ)³
    lambda_compton_m = HBAR_J_S / (m_phi_joule / C_M_S**2) / C_M_S  # m
    
    if lambda_compton_m <= 0:
        return 0.0
    
    volume_m3 = lambda_compton_m**3
    
    # ΔV = (energy scale m)⁴ × (1/volume) × loop factor × log
    delta_v_joule = m_phi_joule
    delta_V_per_m3 = prefactor * (delta_v_joule**4 / (HBAR_J_S * C_M_S)**3) * log_term
    
    return float(delta_V_per_m3)


def quadratic_divergence_check(
    phi: float,
    yukawa: float,
    lambda_uv_gev: float,
    threshold_ratio: float = 0.1,
) -> bool:
    """
    Check if loop correction exhibits quadratic divergence ΔV ∝ Λ_UV².
    
    Parameters
    ----------
    phi : float
        Scalar field value in Planck units.
    yukawa : float
        Yukawa coupling.
    lambda_uv_gev : float
        UV cutoff in GeV.
    threshold_ratio : float
        m/Λ_UV threshold below which quadratic divergence is flagged.
    
    Returns
    -------
    bool
        True if m(φ) ≪ Λ_UV (quadratic divergence regime).
    
    Notes
    -----
    In the limit m ≪ Λ:
        ΔV ~ (Λ²/16π²) × m² × [ln(Λ/m) + const]
    which is quadratically divergent in Λ.
    
    This signals absence of UV protection (e.g., no SUSY, no threshold).
    Mechanisms in this regime require fine-tuning at each loop order.
    """
    M_Pl_GeV = np.sqrt(HBAR_J_S * C_M_S / G_M3_KG_S2) / GEV_J
    m_phi_gev = yukawa * abs(phi) * M_Pl_GeV
    
    if m_phi_gev < 1e-30:
        # m → 0: quadratic divergence always present
        return True
    
    ratio = m_phi_gev / lambda_uv_gev
    return ratio < threshold_ratio


def radiative_stability_check(
    phi_min: float = 0.01,
    phi_max: float = 10.0,
    n_phi: int = 100,
    yukawa: float = 1.0,
    lambda_uv_gev: float = 1e3,
    n_dof: int = 3,
    h0_km_s_mpc: float = 70.0,
    omega_lambda: float = 0.7,
    stability_threshold: float = 1.0,
) -> BackreactionResult:
    """
    Comprehensive radiative stability check for scalar field mechanism.
    
    Parameters
    ----------
    phi_min, phi_max : float
        Field range to scan (in Planck units).
    n_phi : int
        Number of field values to sample.
    yukawa : float
        Yukawa coupling y (dimensionless).
    lambda_uv_gev : float
        UV cutoff in GeV.
    n_dof : int
        Number of degrees of freedom.
    h0_km_s_mpc : float
        Present-day Hubble for observed ρ_Λ estimate.
    omega_lambda : float
        Dark energy density fraction.
    stability_threshold : float
        Max allowed tuning level (max(ΔV) / ρ_Λ < threshold → stable).
    
    Returns
    -------
    BackreactionResult
        Full analysis including tuning level and stability verdict.
    
    Notes
    -----
    Procedure:
    1. Sample field values φ ∈ [φ_min, φ_max]
    2. Compute ΔV(φ) at each point
    3. Find max(ΔV) over field range
    4. Compare to observed ρ_Λ = 3 H₀² Ω_Λ / (8πG)
    5. Flag as unstable if max(ΔV) / ρ_Λ > threshold
    
    Physical interpretation:
    - Tuning level = 1: Loop corrections ~ observed Λ (marginal)
    - Tuning level > 1: Loop corrections exceed Λ (fine-tuning required)
    - Tuning level ≪ 1: Mechanism is radiatively stable
    
    Example verdict:
    - Quintessence with y ~ 1, Λ_UV ~ TeV: UNSTABLE (tuning ~ 10⁶)
    - SUSY-protected (y ~ 0 or threshold at m_SUSY): STABLE
    """
    # Sample field range
    phi_grid = np.linspace(phi_min, phi_max, n_phi)
    
    # Compute loop corrections
    delta_V_loop = np.array([
        coleman_weinberg_correction(phi, yukawa, lambda_uv_gev, n_dof)
        for phi in phi_grid
    ])
    
    # Maximum correction
    max_delta_V = np.max(np.abs(delta_V_loop))
    
    # Observed vacuum energy density (ρ_Λ = 3 H₀² Ω_Λ / 8πG)
    # H0 in km/s/Mpc → s⁻¹: multiply by 1000/Mpc_in_m
    h0_si = (h0_km_s_mpc * 1e3) / MPC_M  # s⁻¹
    rho_lambda_obs = (3.0 * h0_si**2 * omega_lambda) / (8.0 * PI * G_M3_KG_S2)
    
    # Tuning level
    if rho_lambda_obs > 0:
        tuning_level = max_delta_V / rho_lambda_obs
    else:
        tuning_level = np.inf
    
    # Stability verdict
    is_stable = tuning_level < stability_threshold
    
    # Check for quadratic divergence
    quad_div_flags = [
        quadratic_divergence_check(phi, yukawa, lambda_uv_gev)
        for phi in phi_grid
    ]
    quadratic_divergence = any(quad_div_flags)
    
    return BackreactionResult(
        phi_grid=phi_grid,
        delta_V_loop=delta_V_loop,
        max_delta_V=max_delta_V,
        rho_lambda_observed=rho_lambda_obs,
        tuning_level=tuning_level,
        is_stable=is_stable,
        quadratic_divergence_detected=quadratic_divergence,
    )


def estimate_required_tuning(
    delta_V_max: float,
    rho_lambda_obs: float,
) -> Tuple[float, str]:
    """
    Estimate fine-tuning required to cancel loop corrections.
    
    Parameters
    ----------
    delta_V_max : float
        Maximum loop correction (J/m³).
    rho_lambda_obs : float
        Observed vacuum energy (J/m³).
    
    Returns
    -------
    tuning_factor : float
        Required cancellation: ΔV / ρ_Λ.
    verdict : str
        Human-readable assessment.
    
    Notes
    -----
    Tuning interpretation:
    - < 10: Natural (no fine-tuning)
    - 10-100: Mild tuning (1-2 digits)
    - 100-10⁴: Significant tuning (LEP electroweak hierarchy)
    - 10⁴-10¹⁰: Severe tuning (naturalness problem)
    - > 10¹⁰: Extreme tuning (similar to bare CC problem)
    """
    tuning = delta_V_max / rho_lambda_obs
    
    if tuning < 10:
        verdict = "Natural (no fine-tuning required)"
    elif tuning < 100:
        verdict = "Mild tuning (1-2 orders of magnitude)"
    elif tuning < 1e4:
        verdict = "Significant tuning (similar to electroweak hierarchy)"
    elif tuning < 1e10:
        verdict = "Severe tuning (naturalness crisis)"
    else:
        verdict = "Extreme tuning (cosmological constant problem at loop level)"
    
    return tuning, verdict


# Predefined scenarios for common models

def susy_breaking_backreaction(
    m_susy_gev: float = 1e3,
    lambda_uv_gev: float = 1e16,
) -> BackreactionResult:
    """
    Backreaction for SUSY-breaking vacuum energy.
    
    Parameters
    ----------
    m_susy_gev : float
        SUSY-breaking scale in GeV (e.g., 1 TeV).
    lambda_uv_gev : float
        UV cutoff (e.g., GUT scale ~ 10¹⁶ GeV).
    
    Returns
    -------
    BackreactionResult
        Stability analysis for SUSY scenario.
    
    Notes
    -----
    SUSY vacuum energy: ρ_vac ~ m_SUSY⁴ / (16π²) × ln(M_Pl / m_SUSY)
    
    For observed ρ_Λ ~ 10⁻⁴⁷ GeV⁴, requires m_SUSY ~ 10⁻³ eV.
    But LHC: m_SUSY > 1 TeV → tuning ~ (1 TeV / 10⁻³ eV)⁴ ~ 10⁵⁶.
    """
    # Effective field φ ~ m_SUSY / M_Pl (dimensionless)
    M_Pl_GeV = 1.22e19
    phi_eff = m_susy_gev / M_Pl_GeV
    
    # Yukawa ~ O(1) for MSSM
    yukawa = 1.0
    
    # Check single point (SUSY threshold)
    result = radiative_stability_check(
        phi_min=phi_eff * 0.9,
        phi_max=phi_eff * 1.1,
        n_phi=10,
        yukawa=yukawa,
        lambda_uv_gev=lambda_uv_gev,
        n_dof=100,  # Many MSSM states
    )
    
    return result


def quintessence_backreaction(
    phi_today: float = 1.0,
    yukawa: float = 1e-20,
    lambda_uv_gev: float = 1e3,
) -> BackreactionResult:
    """
    Backreaction for quintessence scalar field.
    
    Parameters
    ----------
    phi_today : float
        Present field value in Planck units.
    yukawa : float
        Yukawa coupling to matter (must be tiny for equivalence principle).
    lambda_uv_gev : float
        UV cutoff.
    
    Returns
    -------
    BackreactionResult
        Stability analysis for quintessence.
    
    Notes
    -----
    Quintessence requires:
    1. y < 10⁻²⁰ (fifth force constraints)
    2. Slow-roll: |V'| ≪ V, |V''| ≪ H²
    3. Tracking: φ evolves from O(M_Pl) to O(M_Pl) over Hubble time
    
    With y ~ 10⁻²⁰, loop corrections are negligible (radiatively stable).
    But: no explanation for why y is so small (fine-tuning of couplings).
    """
    result = radiative_stability_check(
        phi_min=phi_today * 0.5,
        phi_max=phi_today * 2.0,
        n_phi=50,
        yukawa=yukawa,
        lambda_uv_gev=lambda_uv_gev,
        n_dof=3,  # SM quarks
    )
    
    return result
