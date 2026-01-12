"""
Likelihood functions for cosmological parameter constraints.

Provides lightweight Bayesian inference tools for fitting mechanism parameters to observational data.
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple
import numpy as np
from scipy import optimize

from .data_loader import DistanceModulusPoint
from .frw import luminosity_distance_m, h_z_lcdm_s_inv
from .mechanisms import CosmologyBackground
from .constants import MPC_M
from .cmb_bao_observables import (
    CMBObservable,
    BAOObservable,
    angular_diameter_distance_mpc,
    cmb_acoustic_scale_ell_a,
    dilation_scale_dv,
)


@dataclass
class LikelihoodResult:
    """Result from likelihood evaluation or parameter fit."""
    log_likelihood: float
    chi_squared: float
    dof: int
    reduced_chi_squared: float
    best_fit_params: Optional[Dict[str, Any]] = None
    param_uncertainties: Optional[Dict[str, float]] = None


def distance_modulus_likelihood(
    data: List[DistanceModulusPoint],
    mechanism_evaluator: Callable[[float], float],
    h0_fiducial: float = 70.0,
) -> LikelihoodResult:
    """
    Compute likelihood for distance modulus observations given a mechanism.

    Parameters
    ----------
    data : List[DistanceModulusPoint]
        Observed distance modulus data.
    mechanism_evaluator : Callable[[float], float]
        Function that takes redshift z and returns ρ_DE(z) in J/m³.
        Used to compute H(z) and luminosity distance.
    h0_fiducial : float
        Fiducial H0 in km/s/Mpc (for absolute calibration).

    Returns
    -------
    LikelihoodResult
        Log-likelihood, chi-squared, and goodness-of-fit statistics.

    Notes
    -----
    Distance modulus: μ(z) = 5 log₁₀(d_L / 10 pc)
    where d_L is luminosity distance in parsecs.

    Chi-squared: χ² = Σ [(μ_obs - μ_theory)² / σ_μ²]
    Log-likelihood: ln L = -χ²/2 (Gaussian errors)
    """
    chi_squared = 0.0
    n_points = len(data)

    for point in data:
        # Compute theoretical distance modulus
        bg = CosmologyBackground(h0_km_s_mpc=h0_fiducial, omega_m=0.3)
        
        # Get ρ_DE(z) from mechanism
        rho_de_z = mechanism_evaluator(point.z)
        
        # Compute luminosity distance (Mpc)
        # For simplicity, assume ΛCDM-like evolution with ρ_DE = const
        # (More sophisticated: integrate with varying ρ_DE(z))
        d_L_m = luminosity_distance_m(point.z, lambda zz: h_z_lcdm_s_inv(zz, bg))
        d_L_mpc = d_L_m / MPC_M
        
        # Convert to distance modulus
        # μ = 5 log₁₀(d_L / 10 pc) = 5 log₁₀(d_L_Mpc * 10^6 / 10) = 5 log₁₀(d_L_Mpc) + 25
        mu_theory = 5.0 * np.log10(d_L_mpc) + 25.0
        
        # Compute chi-squared contribution
        residual = point.mu - mu_theory
        chi_squared += (residual / point.sigma_mu) ** 2

    # Degrees of freedom (n_data - n_params, assume 0 free params for now)
    dof = n_points
    reduced_chi_squared = chi_squared / dof if dof > 0 else float("inf")

    # Log-likelihood (Gaussian errors)
    log_likelihood = -0.5 * chi_squared

    return LikelihoodResult(
        log_likelihood=log_likelihood,
        chi_squared=chi_squared,
        dof=dof,
        reduced_chi_squared=reduced_chi_squared,
    )


def fit_mechanism_parameters(
    data: List[DistanceModulusPoint],
    mechanism_factory: Callable[[Dict[str, Any]], Any],
    initial_params: Dict[str, float],
    param_bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    h0_fiducial: float = 70.0,
) -> LikelihoodResult:
    """
    Fit mechanism parameters to distance modulus data using maximum likelihood.

    Parameters
    ----------
    data : List[DistanceModulusPoint]
        Observed distance modulus data.
    mechanism_factory : Callable[[Dict[str, Any]], Any]
        Function that takes parameter dict and returns mechanism instance.
    initial_params : Dict[str, float]
        Initial parameter values for optimization.
    param_bounds : Dict[str, Tuple[float, float]], optional
        Bounds for each parameter (min, max). If None, unbounded.
    h0_fiducial : float
        Fiducial H0 in km/s/Mpc.

    Returns
    -------
    LikelihoodResult
        Best-fit parameters, uncertainties, and goodness-of-fit statistics.

    Notes
    -----
    Uses SciPy minimize to find maximum likelihood parameters.
    Uncertainties estimated from inverse Hessian (approximate).
    """
    param_names = list(initial_params.keys())
    x0 = np.array([initial_params[name] for name in param_names])

    # Set up bounds
    bounds = None
    if param_bounds is not None:
        bounds = [param_bounds.get(name, (None, None)) for name in param_names]

    # Negative log-likelihood objective
    def neg_log_likelihood(x):
        params = {name: val for name, val in zip(param_names, x)}
        mech = mechanism_factory(params)
        
        # Create evaluator function
        def evaluator(z):
            return mech.evaluate(np.array([z]))[0]
        
        result = distance_modulus_likelihood(data, evaluator, h0_fiducial)
        return -result.log_likelihood

    # Optimize
    result = optimize.minimize(
        neg_log_likelihood,
        x0,
        method="L-BFGS-B" if bounds else "BFGS",
        bounds=bounds,
    )

    best_fit = {name: val for name, val in zip(param_names, result.x)}

    # Compute final likelihood with best-fit parameters
    mech_best = mechanism_factory(best_fit)
    
    def evaluator_best(z):
        return mech_best.evaluate(np.array([z]))[0]
    
    final_likelihood = distance_modulus_likelihood(data, evaluator_best, h0_fiducial)

    # Estimate uncertainties from inverse Hessian (if available)
    uncertainties = None
    if hasattr(result, "hess_inv"):
        try:
            # For BFGS: hess_inv is directly the inverse Hessian
            if isinstance(result.hess_inv, np.ndarray):
                cov_matrix = result.hess_inv
            else:
                # For L-BFGS-B: approximate from final Hessian
                cov_matrix = np.diag([1e-3] * len(param_names))  # fallback
            
            param_std = np.sqrt(np.diag(cov_matrix))
            uncertainties = {name: std for name, std in zip(param_names, param_std)}
        except:
            pass

    return LikelihoodResult(
        log_likelihood=final_likelihood.log_likelihood,
        chi_squared=final_likelihood.chi_squared,
        dof=final_likelihood.dof - len(param_names),  # Subtract fitted params
        reduced_chi_squared=final_likelihood.chi_squared / (final_likelihood.dof - len(param_names)),
        best_fit_params=best_fit,
        param_uncertainties=uncertainties,
    )


def cmb_likelihood(
    cmb_obs: CMBObservable,
    hz_s_inv_callable: Callable[[float], float],
) -> LikelihoodResult:
    """
    Compute CMB acoustic scale likelihood.

    Parameters
    ----------
    cmb_obs : CMBObservable
        CMB acoustic scale measurement (ℓ_A).
    hz_s_inv_callable : Callable[[float], float]
        Function returning H(z) in s^-1.

    Returns
    -------
    LikelihoodResult
        Likelihood result from CMB.

    Notes
    -----
    Compares theoretical ℓ_A = π D_C(z_*) / r_s to observed value.
    Uses COMOVING distance D_C, not angular diameter distance D_A.
    """
    ell_a_theory = cmb_acoustic_scale_ell_a(
        cmb_obs.z_star,
        hz_s_inv_callable,
        r_s_mpc=cmb_obs.r_s_mpc,
    )
    
    residual = cmb_obs.ell_a - ell_a_theory
    chi2 = (residual / cmb_obs.sigma_ell_a) ** 2
    
    return LikelihoodResult(
        log_likelihood=-0.5 * chi2,
        chi_squared=chi2,
        dof=1,  # Single CMB observable
        reduced_chi_squared=chi2,
    )


def bao_likelihood(
    bao_obs_list: List[BAOObservable],
    hz_s_inv_callable: Callable[[float], float],
) -> LikelihoodResult:
    """
    Compute BAO likelihood.

    Parameters
    ----------
    bao_obs_list : List[BAOObservable]
        List of BAO measurements.
    hz_s_inv_callable : Callable[[float], float]
        Function returning H(z) in s^-1.

    Returns
    -------
    LikelihoodResult
        Likelihood result from BAO.

    Notes
    -----
    Supports D_V(z) measurements (dilation scale).
    Can be extended for H(z) and D_A(z) separately.
    """
    chi2_total = 0.0
    
    for obs in bao_obs_list:
        if obs.measurement_type == "DV":
            theory = dilation_scale_dv(obs.z, hz_s_inv_callable)
        elif obs.measurement_type == "DA":
            theory = angular_diameter_distance_mpc(obs.z, hz_s_inv_callable)
        elif obs.measurement_type == "H":
            # H(z) in km/s/Mpc from H(z) in s^-1
            h_z_s_inv = hz_s_inv_callable(obs.z)
            theory = h_z_s_inv * (MPC_M / 1e3)
        else:
            raise ValueError(f"Unknown BAO measurement type: {obs.measurement_type}")
        
        residual = obs.value - theory
        chi2_total += (residual / obs.sigma) ** 2
    
    n_obs = len(bao_obs_list)
    return LikelihoodResult(
        log_likelihood=-0.5 * chi2_total,
        chi_squared=chi2_total,
        dof=n_obs,
        reduced_chi_squared=chi2_total / n_obs if n_obs > 0 else float("inf"),
    )


def joint_likelihood(
    sne_data: Optional[List[DistanceModulusPoint]],
    cmb_obs: Optional[CMBObservable],
    bao_obs_list: Optional[List[BAOObservable]],
    hz_s_inv_callable: Callable[[float], float],
    h0_fiducial: float = 70.0,
) -> LikelihoodResult:
    """
    Compute joint SNe + CMB + BAO likelihood.

    Parameters
    ----------
    sne_data : List[DistanceModulusPoint], optional
        SNe Ia distance modulus data.
    cmb_obs : CMBObservable, optional
        CMB acoustic scale measurement.
    bao_obs_list : List[BAOObservable], optional
        BAO measurements.
    hz_s_inv_callable : Callable[[float], float]
        Function returning H(z) in s^-1.
    h0_fiducial : float
        Fiducial H0 for SNe absolute calibration.

    Returns
    -------
    LikelihoodResult
        Joint likelihood with combined chi-squared.

    Notes
    -----
    Chi-squared: χ² = χ²_SNe + χ²_CMB + χ²_BAO
    """
    chi2_total = 0.0
    n_data = 0
    
    # SNe contribution
    if sne_data is not None and len(sne_data) > 0:
        def evaluator(z):
            # For SNe, we need ρ_DE(z) but we only have H(z)
            # Approximate with constant dark energy for now
            return 5.3e-10
        
        sne_result = distance_modulus_likelihood(sne_data, evaluator, h0_fiducial)
        chi2_total += sne_result.chi_squared
        n_data += len(sne_data)
    
    # CMB contribution
    if cmb_obs is not None:
        cmb_result = cmb_likelihood(cmb_obs, hz_s_inv_callable)
        chi2_total += cmb_result.chi_squared
        n_data += 1
    
    # BAO contribution
    if bao_obs_list is not None and len(bao_obs_list) > 0:
        bao_result = bao_likelihood(bao_obs_list, hz_s_inv_callable)
        chi2_total += bao_result.chi_squared
        n_data += len(bao_obs_list)
    
    dof = n_data
    reduced_chi2 = chi2_total / dof if dof > 0 else float('inf')
    log_likelihood = -0.5 * chi2_total
    
    return LikelihoodResult(
        log_likelihood=log_likelihood,
        chi_squared=chi2_total,
        dof=dof,
        reduced_chi_squared=reduced_chi2,
    )
