# Provenance and scope

This repository is a **workbench**: it encodes common baseline equations and lightweight toy parameterizations for comparison and constraint exercises.

## Baseline equations implemented

- Critical density (mass):
  - $\rho_c = \frac{3H_0^2}{8\pi G}$
- Critical density (energy):
  - $\rho_{c,E} = \rho_c c^2$
- Observed dark-energy density (ΛCDM bookkeeping):
  - $\rho_\Lambda = \Omega_\Lambda\,\rho_{c,E}$
- Cosmological constant parameter conversion:
  - $\Omega_\Lambda = \frac{\Lambda c^2}{3H_0^2}$ so $\Lambda = \frac{3\Omega_\Lambda H_0^2}{c^2}$
- Energy density ↔ Λ conversion:
  - $\rho_\Lambda = \frac{\Lambda c^4}{8\pi G}$

Code locations:
- `ccw.cosmology.observed_lambda_from_h0_omega`
- `ccw.constraints.*`

## Naive QFT scaling (demonstration only)

A common “back-of-envelope” estimate for vacuum energy density with a UV cutoff scales as

- $\rho_{\rm naive} \propto E_{\rm cut}^4$

We use a compact form in code:

- $\rho_{\rm naive} \approx \frac{E_{\rm cut}^4}{16\pi^2\,\hbar^3 c^3}$ for 1 massless DOF.

This is **not** a physical prediction of QFT+GR; it is implemented strictly to reproduce the magnitude mismatch.

Code location:
- `ccw.qft_vacuum.naive_vacuum_energy_density_massless_dof`

## Mechanisms

Mechanism implementations under `ccw.mechanisms.*` are explicitly labeled as **toy/phenomenological** unless stated otherwise. They are designed to be testable, parameter-sweepable, and assumption-traceable.

- CPL: phenomenological w(a) parameterization
- RVM: toy running-vacuum scaling for H(z)
- Unimodular: bookkeeping stub
- Sequestering: residual δρ stub

## What this repo does not claim

- No claim of solving the cosmological constant problem.
- No claim that toy models are derived from a UV-complete theory.
- No claim of observational fits unless/when explicitly added with data provenance.
