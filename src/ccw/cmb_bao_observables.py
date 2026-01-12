"""
Extended observables for CMB and BAO constraints.

Provides angular diameter distance d_A(z) for CMB acoustic scale and
BAO H(z) measurements for joint cosmological constraints.
"""

from dataclasses import dataclass
from typing import Callable, List, Optional
import numpy as np
from scipy import integrate

from .constants import C_M_S, MPC_M
from .frw import comoving_distance_m


@dataclass
class CMBObservable:
    """CMB acoustic scale observable."""
    z_star: float  # Redshift of last scattering (typically ~1090)
    ell_a: float  # Acoustic scale ℓ_A = π D_A / r_s (multipole moment)
    sigma_ell_a: float  # Uncertainty in ℓ_A
    r_s_mpc: float = 147.09  # Sound horizon at decoupling (Mpc, Planck 2018)


@dataclass
class BAOObservable:
    """BAO Hubble parameter or distance measurement."""
    z: float
    measurement_type: str  # "H" or "DV" or "DA"
    value: float  # H(z) in km/s/Mpc, or D_V(z) or D_A(z) in Mpc
    sigma: float  # Uncertainty
    source: str = "BOSS/DESI"


def angular_diameter_distance_mpc(
    z: float,
    hz_s_inv: Callable[[float], float],
) -> float:
    """
    Compute angular diameter distance in Mpc.

    Parameters
    ----------
    z : float
        Redshift.
    hz_s_inv : Callable[[float], float]
        Function returning H(z) in s^-1.

    Returns
    -------
    float
        Angular diameter distance in Mpc.

    Notes
    -----
    For flat FRW: D_A(z) = D_C(z) / (1+z) = comoving_distance / (1+z)
    """
    d_c_m = comoving_distance_m(z, hz_s_inv)
    d_a_m = d_c_m / (1.0 + z)
    return d_a_m / MPC_M


def cmb_acoustic_scale_ell_a(
    z_star: float,
    hz_s_inv: Callable[[float], float],
    r_s_mpc: float = 147.09,
) -> float:
    """
    Compute CMB acoustic scale ℓ_A = π D_C(z_*) / r_s.

    Parameters
    ----------
    z_star : float
        Redshift of last scattering (~1090).
    hz_s_inv : Callable[[float], float]
        Function returning H(z) in s^-1.
    r_s_mpc : float
        Sound horizon at decoupling in Mpc (default: Planck 2018 value).

    Returns
    -------
    float
        Multipole acoustic scale ℓ_A.

    Notes
    -----
    Planck 2018: ℓ_A = 301.63 ± 0.15 (from D_C and r_s).
    Related to first acoustic peak location in CMB power spectrum.
    Uses COMOVING distance D_C, not angular diameter distance D_A.
    """
    # Use comoving distance, not angular diameter distance
    d_c_m = comoving_distance_m(z_star, hz_s_inv)
    d_c_mpc = d_c_m / MPC_M
    ell_a = np.pi * d_c_mpc / r_s_mpc
    return ell_a


def dilation_scale_dv(
    z: float,
    hz_s_inv: Callable[[float], float],
) -> float:
    """
    Compute BAO dilation scale D_V(z) = [(1+z)^2 D_A^2 c z / H(z)]^{1/3}.

    Parameters
    ----------
    z : float
        Redshift.
    hz_s_inv : Callable[[float], float]
        Function returning H(z) in s^-1.

    Returns
    -------
    float
        Dilation scale in Mpc.

    Notes
    -----
    D_V combines line-of-sight and transverse information from BAO.
    Often reported in BAO measurements (e.g., BOSS, DESI).
    """
    d_a_mpc = angular_diameter_distance_mpc(z, hz_s_inv)
    h_z_s_inv = hz_s_inv(z)
    
    # Convert H to km/s/Mpc for dimensional consistency
    h_z_km_s_mpc = h_z_s_inv * (MPC_M / 1e3)
    
    # D_V = [(1+z)^2 D_A^2 c z / H]^{1/3}
    # Units: [Mpc^2 * (km/s) / (km/s/Mpc)]^{1/3} = [Mpc^3]^{1/3} = Mpc
    c_km_s = C_M_S / 1e3
    d_v_cubed = (1 + z)**2 * d_a_mpc**2 * c_km_s * z / h_z_km_s_mpc
    d_v_mpc = d_v_cubed**(1.0/3.0)
    
    return d_v_mpc


def get_planck_cmb_observable() -> CMBObservable:
    """
    Return Planck 2018 CMB acoustic scale measurement.

    Returns
    -------
    CMBObservable
        Planck 2018 ℓ_A measurement.

    Notes
    -----
    Planck Collaboration (2020), "Planck 2018 results. VI. Cosmological parameters"
    ℓ_A = 301.63 ± 0.15 (acoustic scale multipole, Table 2, TT,TE,EE+lowE+lensing)
    z_* = 1089.92 ± 0.25
    r_s(z_*) = (147.09 ± 0.26) Mpc
    """
    return CMBObservable(
        z_star=1089.92,  # Planck 2018
        ell_a=301.63,  # Multipole acoustic scale
        sigma_ell_a=0.15,  # Very precise!
        r_s_mpc=147.09,  # Planck 2018
    )


def get_boss_bao_observables() -> List[BAOObservable]:
    """
    Return BOSS BAO measurements (subset for demonstration).

    Returns
    -------
    List[BAOObservable]
        BOSS DR12 BAO measurements.

    Notes
    -----
    Based on Alam et al. (2017), MNRAS 470, 2617
    Simplified subset for demonstration.
    """
    # BOSS DR12 measurements (simplified)
    measurements = [
        # D_V measurements at effective redshifts
        BAOObservable(z=0.38, measurement_type="DV", value=1512.0, sigma=25.0, source="BOSS DR12"),
        BAOObservable(z=0.51, measurement_type="DV", value=1975.0, sigma=30.0, source="BOSS DR12"),
        BAOObservable(z=0.61, measurement_type="DV", value=2307.0, sigma=40.0, source="BOSS DR12"),
    ]
    
    return measurements


def get_desi_bao_observables() -> List[BAOObservable]:
    """
    Return mock DESI BAO measurements (for demonstration).

    Returns
    -------
    List[BAOObservable]
        Mock DESI BAO measurements.

    Notes
    -----
    DESI will provide percent-level measurements; these are representative values.
    """
    # Mock DESI measurements with improved precision
    measurements = [
        BAOObservable(z=0.30, measurement_type="DV", value=1226.0, sigma=15.0, source="DESI mock"),
        BAOObservable(z=0.50, measurement_type="DV", value=1970.0, sigma=20.0, source="DESI mock"),
        BAOObservable(z=0.70, measurement_type="DV", value=2555.0, sigma=25.0, source="DESI mock"),
        BAOObservable(z=1.00, measurement_type="DV", value=3279.0, sigma=35.0, source="DESI mock"),
    ]
    
    return measurements
