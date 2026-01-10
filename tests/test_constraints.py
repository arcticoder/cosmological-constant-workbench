from __future__ import annotations

import math

from ccw.constraints import lambda_from_rho_lambda, rho_lambda_from_lambda


def test_lambda_rho_roundtrip() -> None:
    rho = 6e-10
    lam = lambda_from_rho_lambda(rho)
    rho2 = rho_lambda_from_lambda(lam)
    assert math.isfinite(lam)
    assert abs(rho2 - rho) / rho < 1e-12
