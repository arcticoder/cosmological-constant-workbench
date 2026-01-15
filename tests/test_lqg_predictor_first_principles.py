"""Tests for configurable external LQG predictor runs.

These are skipped if the optional lqg-cosmological-constant-predictor repo is
not importable.
"""

import math

import pytest

from ccw.integrations.lqg_predictor import (
    LQGPredictorFirstPrinciplesConfig,
    lqg_predictor_available,
    run_first_principles,
)


@pytest.mark.skipif(not lqg_predictor_available(), reason="LQG predictor not available")
def test_run_first_principles_returns_finite_positive_density():
    res = run_first_principles(LQGPredictorFirstPrinciplesConfig(target_scale_m=1e-15, include_uncertainty=False))
    assert math.isfinite(res.lambda_effective_m_minus2)
    assert math.isfinite(res.vacuum_energy_density_j_m3)
    assert res.lambda_effective_m_minus2 > 0
    assert res.vacuum_energy_density_j_m3 > 0


@pytest.mark.skipif(not lqg_predictor_available(), reason="LQG predictor not available")
def test_run_first_principles_accepts_param_overrides():
    # mu_polymer is a documented CLI parameter in the external predictor.
    res = run_first_principles(
        LQGPredictorFirstPrinciplesConfig(
            target_scale_m=1e-15,
            include_uncertainty=False,
            params_overrides={"mu_polymer": 0.2},
        )
    )
    assert res.vacuum_energy_density_j_m3 > 0


@pytest.mark.skipif(not lqg_predictor_available(), reason="LQG predictor not available")
def test_run_first_principles_rejects_unknown_override_keys():
    with pytest.raises(ValueError):
        run_first_principles(
            LQGPredictorFirstPrinciplesConfig(
                target_scale_m=1e-15,
                include_uncertainty=False,
                params_overrides={"definitely_not_a_param": 123},
            )
        )
