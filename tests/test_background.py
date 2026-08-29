"""Acceptance tests for shachen.background (paper section 3.2).

Tolerances: the Planck inverse is analytic, so the roundtrip is exact to
machine precision (rtol=1e-6 leaves headroom). The single golden value was
precomputed independently with scipy.constants (monochromatic Planck at the
10.33 um band center): eps=0.95, T=300 K -> 296.752232 K (rtol=1e-5).
"""

import numpy as np
import pytest
import xarray as xr

from shachen.background import (
    background_bt,
    background_signals,
    planck_radiance,
    planck_temperature,
)
from shachen.constants import BAND_CENTER_UM, Band

TIR_BANDS = (Band.TIR_86, Band.TIR_104, Band.TIR_123)


def test_planck_roundtrip():
    t = np.array([200.0, 240.0, 280.0, 320.0, 340.0])
    for band in TIR_BANDS:
        wl = BAND_CENTER_UM[band]
        np.testing.assert_allclose(planck_temperature(planck_radiance(t, wl), wl), t, rtol=1e-6)


def test_unit_emissivity_is_identity():
    t = np.full((3, 3), 300.0)
    np.testing.assert_allclose(background_bt(t, np.ones((3, 3)), Band.TIR_104), t, atol=1e-6)


def test_background_bt_golden_value():
    bt = background_bt(np.array([300.0]), np.array([0.95]), Band.TIR_104)
    np.testing.assert_allclose(bt, 296.752232, rtol=1e-5)


def test_bt_deficit_grows_with_wavelength():
    # For fixed eps < 1 near 300 K the BT deficit grows with wavelength
    # (toward the Rayleigh-Jeans limit BT -> eps * T).
    deficits = [
        300.0 - float(background_bt(np.array(300.0), np.array(0.90), band)) for band in TIR_BANDS
    ]
    assert deficits[0] < deficits[1] < deficits[2]
    assert all(d > 0 for d in deficits)


def _emissivity(e86, e104, e123, shape=(2, 2)):
    return xr.Dataset(
        {
            "emis_tir_86": xr.DataArray(np.full(shape, e86), dims=("y", "x")),
            "emis_tir_104": xr.DataArray(np.full(shape, e104), dims=("y", "x")),
            "emis_tir_123": xr.DataArray(np.full(shape, e123), dims=("y", "x")),
        }
    )


def _skin(shape=(2, 2), value=300.0):
    return xr.DataArray(np.full(shape, value), dims=("y", "x"))


def test_signals_unit_emissivity_gives_zero_backgrounds():
    out = background_signals(_skin(), _emissivity(1.0, 1.0, 1.0))
    np.testing.assert_allclose(out["rsw_bg"].values, 0.0, atol=1e-6)
    np.testing.assert_allclose(out["btd_bg"].values, 0.0, atol=1e-6)


def test_signals_variable_names():
    out = background_signals(_skin(), _emissivity(1.0, 1.0, 1.0))
    assert {"rsw_bg", "btd_bg", "bt_bg_tir_86", "bt_bg_tir_104", "bt_bg_tir_123"} <= set(
        out.data_vars
    )


def test_signals_nan_emissivity_treated_as_unity():
    # Ocean pixels: CAMEL has no retrieval (NaN) -> eps = 1.0 -> zero background.
    out = background_signals(_skin(), _emissivity(np.nan, np.nan, np.nan))
    np.testing.assert_allclose(out["rsw_bg"].values, 0.0, atol=1e-6)
    np.testing.assert_allclose(out["btd_bg"].values, 0.0, atol=1e-6)
    assert np.isfinite(out["bt_bg_tir_104"].values).all()


def test_signals_desert_like_emissivity():
    # Quartz-like surface: strongly depressed 8.6 um emissivity.
    # Golden values precomputed with the monochromatic Planck formula.
    out = background_signals(_skin(), _emissivity(0.80, 0.96, 0.97))
    np.testing.assert_allclose(out["btd_bg"].values, -8.71202570, rtol=1e-4)
    np.testing.assert_allclose(out["rsw_bg"].values, 0.31103575, rtol=1e-4)
    assert (np.abs(out["btd_bg"].values) > np.abs(out["rsw_bg"].values)).all()


def test_signals_shape_mismatch_raises():
    with pytest.raises(ValueError):
        background_signals(_skin(shape=(3, 3)), _emissivity(1.0, 1.0, 1.0, shape=(2, 2)))
