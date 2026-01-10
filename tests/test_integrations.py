from __future__ import annotations

from ccw.integrations import try_import_coherence_gravity_coupling, try_import_lqg_cosmological_constant_predictor


def test_optional_integrations_are_guarded() -> None:
    # These should never raise; in a multi-repo workspace they may be available.
    r1 = try_import_lqg_cosmological_constant_predictor()
    r2 = try_import_coherence_gravity_coupling()
    assert isinstance(r1.available, bool)
    assert isinstance(r2.available, bool)
