"""Acceptance tests for shachen.dust_tests (Eqs. 13-15).

Pure float64 arithmetic -> rtol=1e-6. Inputs are derived from
DEFAULTS.dust_tests so bound retuning does not invalidate test semantics.

DT3 is asserted in its magnitude-reversed form (see docs/deviations.md, per
the paper's prose): DT3 = clip(((T_MERRA - S) - BT10.4) / depth, 0, 1)
with S applied as printed (T_MERRA - S), S = -10 K land / +5 K ocean.
"""

import numpy as np
import pytest
import xarray as xr

from shachen.constants import DEFAULTS
from shachen.dust_tests import dt1, dt2, dt3, dust_tests

C = DEFAULTS.dust_tests
SHAPE = (2, 2)


def _arr(value, shape=SHAPE, dtype=float):
    return np.full(shape, value, dtype=dtype)


def test_dt1_endpoints_and_midpoint():
    bg = _arr(-1.0)
    np.testing.assert_allclose(dt1(_arr(-1.0), bg), 0.0, atol=1e-12)
    np.testing.assert_allclose(dt1(_arr(C.dt1_max_rsw_k), bg), 1.0, rtol=1e-6)
    mid = (-1.0 + C.dt1_max_rsw_k) / 2.0
    np.testing.assert_allclose(dt1(_arr(mid), bg), 0.5, rtol=1e-6)


def test_dt1_clips():
    bg = _arr(0.0)
    np.testing.assert_allclose(dt1(_arr(C.dt1_max_rsw_k + 2.0), bg), 1.0, rtol=1e-6)
    np.testing.assert_allclose(dt1(_arr(-3.0), bg), 0.0, atol=1e-12)


def test_dt1_degenerate_background_is_zero():
    # No dynamic range left: RSW_bg >= MAX_RSW -> the test cannot fire.
    np.testing.assert_allclose(
        dt1(_arr(C.dt1_max_rsw_k + 2.0), _arr(C.dt1_max_rsw_k)), 0.0, atol=1e-12
    )
    np.testing.assert_allclose(
        dt1(_arr(C.dt1_max_rsw_k + 2.0), _arr(C.dt1_max_rsw_k + 1.0)), 0.0, atol=1e-12
    )


def test_dt2_endpoints_and_midpoint():
    bg = _arr(-0.5)
    np.testing.assert_allclose(dt2(_arr(-0.5), bg), 0.0, atol=1e-12)
    np.testing.assert_allclose(dt2(_arr(C.dt2_max_btd_k), bg), 1.0, rtol=1e-6)
    mid = (-0.5 + C.dt2_max_btd_k) / 2.0
    np.testing.assert_allclose(dt2(_arr(mid), bg), 0.5, rtol=1e-6)
    np.testing.assert_allclose(
        dt2(_arr(C.dt2_max_btd_k + 1.0), _arr(C.dt2_max_btd_k)), 0.0, atol=1e-12
    )


def test_dt3_reversed_land():
    ts = _arr(300.0)
    land = _arr(True, dtype=bool)
    t_ref = 300.0 - C.dt3_shift_land_k  # 310 K with paper defaults
    np.testing.assert_allclose(dt3(_arr(t_ref), ts, land), 0.0, atol=1e-12)
    np.testing.assert_allclose(dt3(_arr(t_ref - C.dt3_depth_k), ts, land), 1.0, rtol=1e-6)
    np.testing.assert_allclose(dt3(_arr(t_ref - C.dt3_depth_k / 2.0), ts, land), 0.5, rtol=1e-6)
    # clips on both sides
    np.testing.assert_allclose(dt3(_arr(t_ref + 10.0), ts, land), 0.0, atol=1e-12)
    np.testing.assert_allclose(dt3(_arr(t_ref - C.dt3_depth_k - 10.0), ts, land), 1.0, rtol=1e-6)


def test_dt3_ocean_shift():
    ts = _arr(300.0)
    ocean = _arr(False, dtype=bool)
    land = _arr(True, dtype=bool)
    t_ref_ocean = 300.0 - C.dt3_shift_ocean_k  # 295 K with paper defaults
    expected_ocean = (t_ref_ocean - 285.0) / C.dt3_depth_k  # 0.2
    np.testing.assert_allclose(dt3(_arr(285.0), ts, ocean), expected_ocean, rtol=1e-6)
    # Same BT over land sits deeper below the (warmer-shifted) reference.
    expected_land = ((300.0 - C.dt3_shift_land_k) - 285.0) / C.dt3_depth_k  # 0.5
    np.testing.assert_allclose(dt3(_arr(285.0), ts, land), expected_land, rtol=1e-6)


def test_dust_tests_assembly():
    scene = xr.Dataset(
        {
            "bt_tir_104": xr.DataArray(_arr(285.0), dims=("y", "x")),
            "bt_tir_123": xr.DataArray(_arr(287.5), dims=("y", "x")),
            "bt_tir_86": xr.DataArray(_arr(286.5), dims=("y", "x")),
        }
    )
    background = xr.Dataset(
        {
            "rsw_bg": xr.DataArray(_arr(0.0), dims=("y", "x")),
            "btd_bg": xr.DataArray(_arr(0.0), dims=("y", "x")),
        }
    )
    ts = xr.DataArray(_arr(300.0), dims=("y", "x"))
    land = xr.DataArray(_arr(True, dtype=bool), dims=("y", "x"))
    out = dust_tests(scene, background, ts, land)
    assert {"dt1", "dt2", "dt3"} <= set(out.data_vars)
    np.testing.assert_allclose(out["dt1"].values, 2.5 / C.dt1_max_rsw_k, rtol=1e-6)
    np.testing.assert_allclose(out["dt2"].values, 1.5 / C.dt2_max_btd_k, rtol=1e-6)
    np.testing.assert_allclose(
        out["dt3"].values, ((300.0 - C.dt3_shift_land_k) - 285.0) / C.dt3_depth_k, rtol=1e-6
    )


def test_nan_propagates():
    obs = _arr(2.0)
    obs[0, 0] = np.nan
    out = dt1(obs, _arr(0.0))
    assert np.isnan(np.asarray(out)[0, 0])
    assert np.isfinite(np.asarray(out)[1, 1])
    bt = _arr(285.0)
    bt[0, 1] = np.nan
    out3 = dt3(bt, _arr(300.0), _arr(True, dtype=bool))
    assert np.isnan(np.asarray(out3)[0, 1])


def test_dust_tests_shape_mismatch_raises():
    scene = xr.Dataset(
        {
            "bt_tir_104": xr.DataArray(_arr(285.0), dims=("y", "x")),
            "bt_tir_123": xr.DataArray(_arr(287.5), dims=("y", "x")),
            "bt_tir_86": xr.DataArray(_arr(286.5), dims=("y", "x")),
        }
    )
    background = xr.Dataset(
        {
            "rsw_bg": xr.DataArray(_arr(0.0), dims=("y", "x")),
            "btd_bg": xr.DataArray(_arr(0.0), dims=("y", "x")),
        }
    )
    ts = xr.DataArray(_arr(300.0), dims=("y", "x"))
    bad_land = xr.DataArray(np.full((3, 3), True), dims=("y", "x"))
    with pytest.raises(ValueError):
        dust_tests(scene, background, ts, bad_land)
