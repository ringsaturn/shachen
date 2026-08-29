"""Acceptance tests for shachen.pipeline.run_debra.

(a) Synthetic end-to-end: a 20x20 all-daytime lat/lon scene over west Texas
    (2017-03-23 18:00 UTC, zenith ~35 deg < 75 -> pure-day blend) with uniform
    ancillary (T_skin = 300 K, eps = 1 -> zero backgrounds, regridding exact),
    split into four quadrants: dust / clear desert / cold cloud / deep
    convection. Every computation is pointwise, so the dust-quadrant CF_comb
    is analytic: DT1 = 2.5/3.5, DT2 = 0.5, DT3 = 0.5, CM_norm_day = 0
    -> CF_comb = N(12/7; 0.25, 2.50) = 41/63 (rtol 1e-6).

(b) Real 2017-03-23 21:45 UTC case (paper Figure 6b proxy), skipped until the
    ABI + MERRA-2 + CAMEL files are all on disk: high mean confidence over the
    W-Texas/NM plume box, near-zero over a non-dust control box.
"""

import datetime as dt
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from shachen.pipeline import run_debra

EXPECTED_VARS = {
    "cf_comb",
    "cf_day",
    "cf_trm",
    "cf_ngt",
    "cm_norm_day",
    "cm_norm_ngt",
    "dt1",
    "dt2",
    "dt3",
    "rsw_bg",
    "btd_bg",
    "zenith_deg",
}

DUST = np.s_[0:10, 0:10]
CLEAR = np.s_[0:10, 10:20]
CLOUD = np.s_[10:20, 0:10]
CONV = np.s_[10:20, 10:20]


def _quadrants(dust, clear, cloud, conv):
    field = np.empty((20, 20))
    field[DUST], field[CLEAR], field[CLOUD], field[CONV] = dust, clear, cloud, conv
    return field


def _synthetic_case(latlon_area):
    scene = xr.Dataset(
        {
            "bt_tir_104": (("y", "x"), _quadrants(285.0, 295.0, 220.0, 210.0)),
            "bt_tir_123": (("y", "x"), _quadrants(287.5, 294.0, 218.0, 209.0)),
            "bt_tir_86": (("y", "x"), _quadrants(286.5, 294.5, 219.0, 209.5)),
            "bt_wv_62": (("y", "x"), _quadrants(245.0, 255.0, 215.0, 209.0)),
            "bt_swir_39": (("y", "x"), _quadrants(286.0, 300.0, 222.0, 211.0)),
        },
        attrs={"area": latlon_area, "start_time": dt.datetime(2017, 3, 23, 18, 0)},
    )
    ts = xr.DataArray(
        np.full((16, 16), 300.0),
        coords={"lat": np.arange(25.0, 41.0), "lon": np.arange(-110.0, -94.0)},
        dims=("lat", "lon"),
    )
    emis = xr.Dataset(
        {
            f"emis_{band}": xr.DataArray(
                np.ones((16, 16)),
                coords={
                    "latitude": np.arange(25.0, 41.0),
                    "longitude": np.arange(-110.0, -94.0),
                },
                dims=("latitude", "longitude"),
            )
            for band in ("tir_86", "tir_104", "tir_123")
        }
    )
    return scene, ts, emis


@pytest.fixture
def synthetic_output(latlon_area):
    scene, ts, emis = _synthetic_case(latlon_area)
    return run_debra(scene, ts, emis)


def test_output_variables_and_attrs(latlon_area, synthetic_output):
    out = synthetic_output
    assert EXPECTED_VARS <= set(out.data_vars)
    for name in EXPECTED_VARS:
        assert out[name].shape == (20, 20)
    assert out.attrs["area"] is latlon_area
    assert out.attrs["start_time"] == dt.datetime(2017, 3, 23, 18, 0)


def test_scene_is_pure_day(synthetic_output):
    # 18:00 UTC over west Texas: zenith ~35 deg everywhere -> pure-day blend.
    assert (synthetic_output["zenith_deg"].values < 75.0).all()


def test_dust_quadrant_analytic(synthetic_output):
    out = synthetic_output
    np.testing.assert_allclose(out["dt1"].values[DUST], 2.5 / 3.5, rtol=1e-6)
    np.testing.assert_allclose(out["dt2"].values[DUST], 0.5, rtol=1e-6)
    np.testing.assert_allclose(out["dt3"].values[DUST], 0.5, rtol=1e-6)
    np.testing.assert_allclose(out["rsw_bg"].values, 0.0, atol=1e-6)
    np.testing.assert_allclose(out["cm_norm_day"].values[DUST], 0.0, atol=1e-12)
    np.testing.assert_allclose(out["cf_comb"].values[DUST], 41.0 / 63.0, rtol=1e-6)


def test_non_dust_quadrants_near_zero(synthetic_output):
    cf = synthetic_output["cf_comb"].values
    assert (cf[CLEAR] < 0.05).all()
    assert (cf[CLOUD] < 0.05).all()
    assert (cf[CONV] < 0.05).all()


def test_cf_comb_in_unit_interval(synthetic_output):
    cf = synthetic_output["cf_comb"].values
    assert np.isfinite(cf).all()
    assert ((cf >= 0.0) & (cf <= 1.0)).all()


# --- precomputed background ---------------------------------------------------


def _uniform_background(shape=(20, 20)):
    """Scene-grid background equivalent to the T=300 K / eps=1 semianalytic one."""
    data = {
        name: (("y", "x"), np.full(shape, 300.0))
        for name in ("bt_bg_tir_86", "bt_bg_tir_104", "bt_bg_tir_123")
    }
    data["rsw_bg"] = (("y", "x"), np.zeros(shape))
    data["btd_bg"] = (("y", "x"), np.zeros(shape))
    data["n_valid"] = (("y", "x"), np.full(shape, 14, dtype=np.int16))
    return xr.Dataset(data)


def test_precomputed_background_matches_semianalytic(latlon_area, synthetic_output):
    scene, ts, _ = _synthetic_case(latlon_area)
    out = run_debra(scene, ts, background=_uniform_background())
    np.testing.assert_allclose(out["cf_comb"].values, synthetic_output["cf_comb"].values, rtol=1e-6)
    np.testing.assert_allclose(out["rsw_bg"].values, 0.0, atol=1e-12)
    assert "n_valid" in out.data_vars
    assert (out["n_valid"].values == 14).all()


def test_exactly_one_background_source_required(latlon_area):
    scene, ts, emis = _synthetic_case(latlon_area)
    with pytest.raises(ValueError):
        run_debra(scene, ts)  # neither source
    with pytest.raises(ValueError):
        run_debra(scene, ts, emis, background=_uniform_background())  # both


def test_background_shape_mismatch_raises(latlon_area):
    scene, ts, _ = _synthetic_case(latlon_area)
    with pytest.raises(ValueError):
        run_debra(scene, ts, background=_uniform_background(shape=(10, 10)))


def test_background_missing_variable_raises(latlon_area):
    scene, ts, _ = _synthetic_case(latlon_area)
    with pytest.raises(ValueError):
        run_debra(scene, ts, background=_uniform_background().drop_vars("rsw_bg"))


# --- (b) real case: paper Figure 6b proxy -----------------------------------

_DATA = Path(__file__).resolve().parents[1] / "data"
_ABI = _DATA / "2017-03-23-swus" / "abi"
_MERRA = _DATA / "merra2" / "merra2_ts_20170323.nc"
_EMIS = _DATA / "emissivity" / "CAM5K30EM_201703.nc"
_HAVE_CASE = (
    _ABI.exists() and len(list(_ABI.glob("*.nc"))) >= 7 and _MERRA.exists() and _EMIS.exists()
)


@pytest.mark.skipif(
    not _HAVE_CASE,
    reason="2017-03-23 case incomplete (need ABI + MERRA-2 + CAMEL; "
    "run scripts/fetch_case.py after the Earthdata password reset)",
)
def test_real_case_figure6b_proxy():
    from shachen.io.emissivity import load_band_emissivity
    from shachen.io.merra import load_skin_temperature
    from shachen.io.satellite import load_scene

    scene = load_scene(sorted(_ABI.glob("*.nc")))
    ts = load_skin_temperature(_MERRA, scene.attrs["start_time"])
    emis = load_band_emissivity(_EMIS)
    out = run_debra(scene, ts, emis)

    cf = np.asarray(out["cf_comb"].values, dtype=float)
    assert cf.shape == scene["bt_tir_104"].shape
    finite = np.isfinite(cf)
    assert finite.mean() > 0.5
    assert ((cf[finite] >= 0.0) & (cf[finite] <= 1.0)).all()

    lons, lats = scene.attrs["area"].get_lonlats()
    # Plume core at 21:45 UTC sits over the El Paso / south-central New Mexico
    # corridor (verified from the scene's own positive-RSW field: peak cells at
    # 30-33N, 104-107W with mean RSW up to +3.9 K), not the W-Texas plains.
    plume = finite & (lats >= 30.0) & (lats <= 33.0) & (lons >= -107.0) & (lons <= -104.5)
    control = finite & (lats >= 24.0) & (lats <= 27.0) & (lons >= -90.0) & (lons <= -86.0)
    assert plume.any() and control.any()
    plume_mean = cf[plume].mean()
    control_mean = cf[control].mean()
    assert plume_mean > 0.3  # confident dust over the W-Texas/NM plume (Fig. 6b)
    assert plume_mean > 3.0 * control_mean
