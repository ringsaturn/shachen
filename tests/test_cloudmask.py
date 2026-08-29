"""Acceptance tests for shachen.cloudmask (Eqs. 1-12).

All expected values follow from the equations evaluated at endpoints and
midpoints of the DEFAULTS bounds (pure float64 arithmetic -> rtol=1e-6).
Test inputs are derived from DEFAULTS.cloud_mask so a retuning of the bounds
does not invalidate the test semantics.

Two deliberate deviations from the printed equations (see docs/deviations.md and
the module docstring): CM2 is magnitude-reversed (Eq. 4 as printed would
saturate the mask in clear sky), and CM_day uses CM3, not the misprinted CM4.
"""

import numpy as np
import pytest
import xarray as xr

from shachen.cloudmask import cloud_mask
from shachen.constants import DEFAULTS

C = DEFAULTS.cloud_mask
TSKIN = 300.0
SHAPE = (4, 4)

EXPECTED_VARS = {
    "cm1",
    "cm2",
    "cm3",
    "cm4",
    "r1",
    "r2_day",
    "r2_ngt",
    "cm_day",
    "cm_ngt",
    "cm_norm_day",
    "cm_norm_ngt",
}


def _mid(bounds):
    return (bounds.min + bounds.max) / 2.0


def _scene(bt104, bt123=None, bt86=None, bt62=None, bt39=None):
    """Uniform scene; unset bands default to values that zero the other tests.

    bt123 = bt104 -> CM3 = 0 and R1's N-term = 0; bt86 = bt104 + r2.min -> R2's
    N-term = 0; bt62 = bt104 - (cm2.max + 15) -> CM2 = 0 (reversed form);
    bt39 = bt104 + cm4.min -> CM4 = 0.
    """
    if bt123 is None:
        bt123 = bt104
    if bt86 is None:
        bt86 = bt104 + C.r2.min
    if bt62 is None:
        bt62 = bt104 - (C.cm2.max + 15.0)
    if bt39 is None:
        bt39 = bt104 + C.cm4.min
    data = {
        "bt_tir_104": bt104,
        "bt_tir_123": bt123,
        "bt_tir_86": bt86,
        "bt_wv_62": bt62,
        "bt_swir_39": bt39,
    }
    return xr.Dataset(
        {k: xr.DataArray(np.full(SHAPE, v, dtype=float), dims=("y", "x")) for k, v in data.items()}
    )


def _skin(value=TSKIN, shape=SHAPE):
    return xr.DataArray(np.full(shape, value), dims=("y", "x"))


def test_output_variables_present():
    out = cloud_mask(_scene(bt104=TSKIN), _skin())
    assert EXPECTED_VARS <= set(out.data_vars)


def test_cm1_endpoints_and_midpoint():
    # Eq. 1: 1 - N(BT10.4; T_skin - offset, T_skin)
    warm = cloud_mask(_scene(bt104=TSKIN), _skin())
    np.testing.assert_allclose(warm["cm1"].values, 0.0, atol=1e-12)
    cold = cloud_mask(_scene(bt104=TSKIN - C.cm1_cold_offset_k), _skin())
    np.testing.assert_allclose(cold["cm1"].values, 1.0, rtol=1e-6)
    half = cloud_mask(_scene(bt104=TSKIN - C.cm1_cold_offset_k / 2.0), _skin())
    np.testing.assert_allclose(half["cm1"].values, 0.5, rtol=1e-6)


def test_cm2_is_magnitude_reversed():
    # Deliberate deviation: CM2 = 1 - N(BT10.4 - BT6.2; cm2 bounds).
    bt104 = TSKIN
    clear = cloud_mask(_scene(bt104, bt62=bt104 - (C.cm2.max + 15.0)), _skin())
    np.testing.assert_allclose(clear["cm2"].values, 0.0, atol=1e-12)
    conv = cloud_mask(_scene(bt104, bt62=bt104 - C.cm2.min), _skin())
    np.testing.assert_allclose(conv["cm2"].values, 1.0, rtol=1e-6)
    half = cloud_mask(_scene(bt104, bt62=bt104 - _mid(C.cm2)), _skin())
    np.testing.assert_allclose(half["cm2"].values, 0.5, rtol=1e-6)


def test_cm3_split_window():
    bt104 = TSKIN
    lo = cloud_mask(_scene(bt104, bt123=bt104 - C.cm3.min), _skin())
    np.testing.assert_allclose(lo["cm3"].values, 0.0, atol=1e-12)
    hi = cloud_mask(_scene(bt104, bt123=bt104 - C.cm3.max), _skin())
    np.testing.assert_allclose(hi["cm3"].values, 1.0, rtol=1e-6)
    half = cloud_mask(_scene(bt104, bt123=bt104 - _mid(C.cm3)), _skin())
    np.testing.assert_allclose(half["cm3"].values, 0.5, rtol=1e-6)


def test_cm4_night_cirrus():
    bt104 = TSKIN
    half = cloud_mask(_scene(bt104, bt39=bt104 + _mid(C.cm4)), _skin())
    np.testing.assert_allclose(half["cm4"].values, 0.5, rtol=1e-6)
    hi = cloud_mask(_scene(bt104, bt39=bt104 + C.cm4.max), _skin())
    np.testing.assert_allclose(hi["cm4"].values, 1.0, rtol=1e-6)


def test_r1_restoral_and_cm1_weighting():
    # Eq. 7 (erratum): R1 = N(BT12.3 - BT10.4; r1) * (1 - CM1)
    warm = cloud_mask(_scene(TSKIN, bt123=TSKIN + _mid(C.r1)), _skin())
    np.testing.assert_allclose(warm["r1"].values, 0.5, rtol=1e-6)
    cold = cloud_mask(
        _scene(TSKIN - C.cm1_cold_offset_k, bt123=TSKIN - C.cm1_cold_offset_k + _mid(C.r1)),
        _skin(),
    )
    np.testing.assert_allclose(cold["r1"].values, 0.0, atol=1e-12)


def test_r2_restorals_and_weights():
    bt104 = TSKIN
    # Eq. 8: N-term at midpoint = 0.5, CM2 = CM3 = 0 -> r2_day = 0.5
    base = cloud_mask(_scene(bt104, bt86=bt104 + _mid(C.r2)), _skin())
    np.testing.assert_allclose(base["r2_day"].values, 0.5, rtol=1e-6)
    # Deep convection (CM2 = 1) removes the day restoral
    conv = cloud_mask(_scene(bt104, bt86=bt104 + _mid(C.r2), bt62=bt104 - C.cm2.min), _skin())
    np.testing.assert_allclose(conv["r2_day"].values, 0.0, atol=1e-12)
    # Eq. 9: CM4 = 1 removes the night restoral but not the day one
    cirrus = cloud_mask(_scene(bt104, bt86=bt104 + _mid(C.r2), bt39=bt104 + C.cm4.max), _skin())
    np.testing.assert_allclose(cirrus["r2_ngt"].values, 0.0, atol=1e-12)
    np.testing.assert_allclose(cirrus["r2_day"].values, 0.5, rtol=1e-6)


def test_day_mask_uses_cm3_night_mask_uses_cm4():
    # Deliberate deviation: Eq. 11 uses CM3 (printed CM4 is a typo).
    bt104 = TSKIN
    only_cm4 = cloud_mask(_scene(bt104, bt39=bt104 + C.cm4.max), _skin())
    np.testing.assert_allclose(only_cm4["cm_day"].values, 0.0, atol=1e-12)
    np.testing.assert_allclose(only_cm4["cm_ngt"].values, 1.0, rtol=1e-6)
    only_cm3 = cloud_mask(_scene(bt104, bt123=bt104 - C.cm3.max), _skin())
    np.testing.assert_allclose(only_cm3["cm_day"].values, 1.0, rtol=1e-6)
    np.testing.assert_allclose(only_cm3["cm_ngt"].values, 0.0, atol=1e-12)


def test_composed_dusty_scene_golden():
    # Dust-like pixel: warm-ish, positive split window, restorals active.
    # Hand-derived: CM1 = 0.3, CM2 = CM3 = 0, R1 = N(2.5; 0, 3.5)*0.7 = 0.5,
    # R2_day = N(1.5; -1, 3) = 0.625 -> CM_day = 0.3 * 0.375 = 0.1125 -> norm 0.
    out = cloud_mask(_scene(285.0, bt123=287.5, bt86=286.5, bt62=245.0, bt39=286.0), _skin())
    np.testing.assert_allclose(out["cm1"].values, 0.3, rtol=1e-6)
    np.testing.assert_allclose(out["r1"].values, 0.5, rtol=1e-6)
    np.testing.assert_allclose(out["r2_day"].values, 0.625, rtol=1e-6)
    np.testing.assert_allclose(out["cm_day"].values, 0.1125, rtol=1e-6)
    np.testing.assert_allclose(out["cm_norm_day"].values, 0.0, atol=1e-12)


def test_cm_norm_midpoint():
    # CM1 alone = midpoint of the cm_norm bounds -> CM_norm = 0.5 (Eq. 12).
    bt104 = TSKIN - C.cm1_cold_offset_k * _mid(C.cm_norm)
    out = cloud_mask(_scene(bt104), _skin())
    np.testing.assert_allclose(out["cm_day"].values, _mid(C.cm_norm), rtol=1e-6)
    np.testing.assert_allclose(out["cm_norm_day"].values, 0.5, rtol=1e-6)


def test_nan_propagates():
    scene = _scene(bt104=TSKIN)
    bt = scene["bt_tir_104"].values.copy()
    bt[1, 2] = np.nan
    scene["bt_tir_104"] = xr.DataArray(bt, dims=("y", "x"))
    out = cloud_mask(scene, _skin())
    assert np.isnan(out["cm1"].values[1, 2])
    assert np.isnan(out["cm_norm_day"].values[1, 2])
    assert np.isfinite(out["cm_norm_day"].values[0, 0])


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        cloud_mask(_scene(bt104=TSKIN), _skin(shape=(3, 3)))
