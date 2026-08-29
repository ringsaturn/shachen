"""Acceptance tests for debra.confidence (Eqs. 16-22, erratum forms of 21-22).

Pure float64 arithmetic -> rtol=1e-6, except assertions involving the
trigonometric blend weights, which use precomputed goldens at rtol=1e-4:
b_ngt_trm(97.5 deg) = ((cos 97.5 - cos 105) / (cos 90 - cos 105))^1.5
              = 0.348987 (NOAA-independent hand computation).
Inputs derive from DEFAULTS.confidence where bounds are involved.
"""

import numpy as np
import pytest
import xarray as xr

from shachen.confidence import confidence
from shachen.constants import DEFAULTS

C = DEFAULTS.confidence
SHAPE = (2, 2)

B_NGT_TRM_97_5 = 0.348987

EXPECTED_VARS = {"cf_day", "cf_trm", "cf_ngt", "b_ngt_trm", "b_trm_day", "cf_comb"}


def _inputs(d1, d2, d3, cm_day=0.0, cm_ngt=0.0, zen=30.0):
    tests = xr.Dataset(
        {
            "dt1": xr.DataArray(np.full(SHAPE, d1), dims=("y", "x")),
            "dt2": xr.DataArray(np.full(SHAPE, d2), dims=("y", "x")),
            "dt3": xr.DataArray(np.full(SHAPE, d3), dims=("y", "x")),
        }
    )
    cloud = xr.Dataset(
        {
            "cm_norm_day": xr.DataArray(np.full(SHAPE, cm_day), dims=("y", "x")),
            "cm_norm_ngt": xr.DataArray(np.full(SHAPE, cm_ngt), dims=("y", "x")),
        }
    )
    zenith = xr.DataArray(np.full(SHAPE, zen), dims=("y", "x"))
    return tests, cloud, zenith


def _norm(raw):
    return np.clip((raw - C.cf_norm.min) / (C.cf_norm.max - C.cf_norm.min), 0.0, 1.0)


def test_output_variables_present():
    out = confidence(*_inputs(0.5, 0.5, 0.5))
    assert EXPECTED_VARS <= set(out.data_vars)


def test_full_day_saturates():
    # Eq. 16 with all tests = 1, no cloud: raw = 3.0 -> Eq. 19 norm -> 1.0
    out = confidence(*_inputs(1.0, 1.0, 1.0, cm_day=0.0, zen=30.0))
    np.testing.assert_allclose(out["b_trm_day"].values, 1.0, rtol=1e-6)
    np.testing.assert_allclose(out["cf_day"].values, 1.0, rtol=1e-6)
    np.testing.assert_allclose(out["cf_comb"].values, 1.0, rtol=1e-6)


def test_zero_signal_is_zero():
    out = confidence(*_inputs(0.0, 0.0, 0.0, zen=30.0))
    np.testing.assert_allclose(out["cf_comb"].values, 0.0, atol=1e-12)


def test_day_norm_midpoint():
    # DT sum at the midpoint of the cf_norm bounds -> cf_day = 0.5 exactly.
    mid = (C.cf_norm.min + C.cf_norm.max) / 2.0
    out = confidence(*_inputs(mid - 1.0, 0.5, 0.5, cm_day=0.0, zen=30.0))
    np.testing.assert_allclose(out["cf_comb"].values, 0.5, rtol=1e-6)


def test_night_uses_max_and_night_mask():
    # Eq. 18: max(DT1, DT2) + 0.5*DT3, suppressed by cm_norm_ngt (not _day).
    out = confidence(*_inputs(0.2, 0.8, 0.4, cm_day=1.0, cm_ngt=0.0, zen=120.0))
    raw = max(0.2, 0.8) + C.dt3_weight_ngt * 0.4  # = 1.0
    np.testing.assert_allclose(out["b_trm_day"].values, 0.0, atol=1e-12)
    np.testing.assert_allclose(out["b_ngt_trm"].values, 0.0, atol=1e-12)
    np.testing.assert_allclose(out["cf_comb"].values, _norm(raw), rtol=1e-6)


def test_terminator_selects_cf_trm():
    # At zenith 90: b_trm_day = 0, b_ngt_trm = 1 -> cf_comb == cf_trm (Eq. 17).
    out = confidence(*_inputs(0.6, 0.2, 0.8, cm_day=0.0, zen=90.0))
    raw = 0.6 + 0.2 + C.dt3_weight_trm * 0.8  # = 1.2
    np.testing.assert_allclose(out["cf_trm"].values, _norm(raw), rtol=1e-6)
    np.testing.assert_allclose(out["cf_comb"].values, out["cf_trm"].values, rtol=1e-6)


def test_blend_weight_golden_and_composition():
    # Eq. 22 (erratum) at zenith 97.5 deg with distinct cf_day/trm/ngt.
    out = confidence(*_inputs(0.6, 0.2, 0.8, cm_day=0.0, cm_ngt=0.5, zen=97.5))
    np.testing.assert_allclose(out["b_ngt_trm"].values, B_NGT_TRM_97_5, rtol=1e-4)
    np.testing.assert_allclose(out["b_trm_day"].values, 0.0, atol=1e-12)
    cf_trm = _norm(0.6 + 0.2 + C.dt3_weight_trm * 0.8)
    cf_ngt = _norm((max(0.6, 0.2) + C.dt3_weight_ngt * 0.8) * (1.0 - 0.5))
    expected = B_NGT_TRM_97_5 * cf_trm + (1.0 - B_NGT_TRM_97_5) * cf_ngt
    np.testing.assert_allclose(out["cf_comb"].values, expected, rtol=1e-4)


def test_cloud_suppression():
    out = confidence(*_inputs(1.0, 1.0, 1.0, cm_day=1.0, zen=30.0))
    np.testing.assert_allclose(out["cf_comb"].values, 0.0, atol=1e-12)


def test_cf_comb_in_unit_interval():
    for zen in (30.0, 80.0, 90.0, 100.0, 120.0):
        out = confidence(*_inputs(0.9, 0.7, 0.6, cm_day=0.2, cm_ngt=0.1, zen=zen))
        v = out["cf_comb"].values
        assert ((v >= 0.0) & (v <= 1.0)).all()


def test_nan_propagates():
    tests, cloud, zenith = _inputs(0.5, 0.5, 0.5)
    z = zenith.values.copy()
    z[0, 0] = np.nan
    out = confidence(tests, cloud, xr.DataArray(z, dims=("y", "x")))
    assert np.isnan(out["cf_comb"].values[0, 0])
    assert np.isfinite(out["cf_comb"].values[1, 1])
    d = tests["dt3"].values.copy()
    d[1, 0] = np.nan
    tests["dt3"] = xr.DataArray(d, dims=("y", "x"))
    out2 = confidence(tests, cloud, zenith)
    assert np.isnan(out2["cf_comb"].values[1, 0])


def test_shape_mismatch_raises():
    tests, cloud, _ = _inputs(0.5, 0.5, 0.5)
    bad_zen = xr.DataArray(np.full((3, 3), 30.0), dims=("y", "x"))
    with pytest.raises(ValueError):
        confidence(tests, cloud, bad_zen)
