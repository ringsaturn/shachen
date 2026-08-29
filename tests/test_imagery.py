"""Acceptance tests for shachen.imagery (Eqs. 23-29, erratum Eqs. 24-25).

All fields are synthetic and every expectation is an exact algebraic value of
the printed equations, so tolerances are tight (rtol 1e-9); the design doc
defines no looser physical tolerance.
"""

import numpy as np
import pytest
import xarray as xr

from shachen.constants import COLOR_DIMMING
from shachen.imagery import baseline_image, debra_imagery, enhanced_rgb, to_uint8

RTOL = 1e-9


def _da(values) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=float), dims=("y", "x"))


def _scene(refl, bt) -> xr.Dataset:
    return xr.Dataset({"refl_vis_064": _da(refl), "bt_tir_104": _da(bt)})


# --- baseline image: Eqs. 23-26 ----------------------------------------------


class TestBaselineImage:
    def test_pure_day_is_domain_normalized_vis(self):
        # Eq. 23 + Eq. 26 with b_bg = 1 (zenith 30 deg << 79 deg).
        scene = _scene([[0.0, 10.0], [20.0, 30.0]], [[250.0, 260.0], [280.0, 300.0]])
        out = baseline_image(scene, _da(np.full((2, 2), 30.0)))
        expected = [[0.0, 1.0 / 3.0], [2.0 / 3.0, 1.0]]
        np.testing.assert_allclose(out["vis_bg"].values, expected, rtol=RTOL)
        np.testing.assert_allclose(out["b_bg"].values, 1.0, rtol=RTOL)
        np.testing.assert_allclose(out["bi"].values, expected, rtol=RTOL)

    def test_pure_night_is_inverted_ir(self):
        # Eq. 24 (erratum): 1 - N(BT; Tmin, Tmax); zenith 120 deg -> b_bg = 0.
        scene = _scene(np.zeros((2, 2)), [[250.0, 260.0], [280.0, 300.0]])
        out = baseline_image(scene, _da(np.full((2, 2), 120.0)))
        expected = [[1.0, 0.8], [0.4, 0.0]]
        np.testing.assert_allclose(out["ir_bg"].values, expected, rtol=RTOL)
        np.testing.assert_allclose(out["b_bg"].values, 0.0, atol=1e-12)
        np.testing.assert_allclose(out["bi"].values, expected, rtol=RTOL)

    def test_b_bg_ramp(self):
        # Eq. 25 (erratum): 1 - N(zenith; 79, 89)**1.5, zenith-degree space.
        zen = [[70.0, 79.0, 84.0, 89.0, 100.0]]
        scene = _scene([[0.0, 1.0, 2.0, 3.0, 4.0]], [[250.0, 260.0, 270.0, 280.0, 290.0]])
        out = baseline_image(scene, _da(zen))
        expected = [[1.0, 1.0, 1.0 - 0.5**1.5, 0.0, 0.0]]
        np.testing.assert_allclose(out["b_bg"].values, expected, rtol=RTOL, atol=1e-12)

    def test_zero_weight_component_ignores_nan(self):
        # Eq. 26: a NaN component with exactly zero blend weight must not
        # poison bi; with nonzero weight NaN propagates.
        zen = [[30.0, 120.0, 85.0], [30.0, 30.0, 30.0]]
        refl = [[0.0, np.nan, np.nan], [np.nan, 1.0, 0.5]]
        bt = [[250.0, 260.0, 270.0], [280.0, 290.0, 300.0]]
        out = baseline_image(_scene(refl, bt), _da(zen))
        bi = out["bi"].values
        np.testing.assert_allclose(bi[0, 0], 0.0, atol=1e-12)  # day, vis valid
        np.testing.assert_allclose(bi[0, 1], 0.8, rtol=RTOL)  # night: ir_bg only
        assert np.isnan(bi[0, 2])  # terminator: vis weight > 0
        assert np.isnan(bi[1, 0])  # day: vis weight > 0
        np.testing.assert_allclose(bi[1, 1], 1.0, rtol=RTOL)
        np.testing.assert_allclose(bi[1, 2], 0.5, rtol=RTOL)

    def test_degenerate_vis_domain_is_zero(self):
        # Constant (or all-NaN) VIS: component defined as 0.0, no error.
        scene = _scene(np.zeros((2, 2)), [[250.0, 260.0], [280.0, 300.0]])
        out = baseline_image(scene, _da(np.full((2, 2), 30.0)))
        np.testing.assert_allclose(out["vis_bg"].values, 0.0, atol=1e-12)
        np.testing.assert_allclose(out["bi"].values, 0.0, atol=1e-12)

        scene_nan = _scene(np.full((2, 2), np.nan), [[250.0, 260.0], [280.0, 300.0]])
        out_nan = baseline_image(scene_nan, _da(np.full((2, 2), 120.0)))
        np.testing.assert_allclose(out_nan["vis_bg"].values, 0.0, atol=1e-12)
        np.testing.assert_allclose(out_nan["bi"].values, [[1.0, 0.8], [0.4, 0.0]], rtol=RTOL)

    def test_bi_in_unit_interval(self):
        rng_refl = [[5.0, 55.0], [25.0, 90.0]]
        bt = [[250.0, 300.0], [270.0, 285.0]]
        out = baseline_image(_scene(rng_refl, bt), _da([[30.0, 85.0], [95.0, 120.0]]))
        bi = out["bi"].values
        assert np.isfinite(bi).all()
        assert ((bi >= 0.0) & (bi <= 1.0)).all()

    def test_missing_variable_raises(self):
        zen = _da(np.full((2, 2), 30.0))
        vis_only = xr.Dataset({"refl_vis_064": _da(np.zeros((2, 2)))})
        ir_only = xr.Dataset({"bt_tir_104": _da(np.full((2, 2), 280.0))})
        with pytest.raises(ValueError):
            baseline_image(vis_only, zen)
        with pytest.raises(ValueError):
            baseline_image(ir_only, zen)

    def test_shape_mismatch_raises(self):
        scene = _scene(np.zeros((2, 2)), np.full((2, 2), 280.0))
        with pytest.raises(ValueError):
            baseline_image(scene, _da(np.full((3, 3), 30.0)))


# --- color modulation: Eqs. 27-29 + gun rescale ------------------------------


class TestEnhancedRgb:
    def test_dims_and_coords(self):
        rgb = enhanced_rgb(_da([[0.5]]), _da([[0.0]]))
        assert rgb.dims == ("y", "x", "gun")
        assert list(np.asarray(rgb.coords["gun"].values)) == ["r", "g", "b"]

    def test_cf_zero_is_grayscale(self):
        # Eqs. 27-29 with CF = 0: every gun is N(BI; 0, 1.2) = BI/1.2.
        bi = _da([[0.0, 0.6], [1.0, 0.3]])
        rgb = enhanced_rgb(bi, _da(np.zeros((2, 2))))
        for gun in range(3):
            np.testing.assert_allclose(rgb.values[..., gun], bi.values / 1.2, rtol=RTOL, atol=1e-12)

    def test_yellow_saturated_dust(self):
        # BI = 0.5, CF = 1.0: red = grn = (0.5*0.5 + 1)/1.2 -> clipped to 1.0;
        # blu = (0.25 + 0.10)/1.2.
        rgb = enhanced_rgb(_da([[0.5]]), _da([[1.0]]))
        np.testing.assert_allclose(rgb.values[0, 0, 0], 1.0, rtol=RTOL)
        np.testing.assert_allclose(rgb.values[0, 0, 1], 1.0, rtol=RTOL)
        np.testing.assert_allclose(rgb.values[0, 0, 2], 0.35 / 1.2, rtol=RTOL)

    def test_cf_below_cap(self):
        # BI = 0.8, CF = 0.4 < cf_cap: red = grn = (0.8*0.6 + 0.4)/1.2,
        # blu = (0.48 + 0.04)/1.2.
        rgb = enhanced_rgb(_da([[0.8]]), _da([[0.4]]))
        np.testing.assert_allclose(rgb.values[0, 0, 0], 0.88 / 1.2, rtol=RTOL)
        np.testing.assert_allclose(rgb.values[0, 0, 1], 0.88 / 1.2, rtol=RTOL)
        np.testing.assert_allclose(rgb.values[0, 0, 2], 0.52 / 1.2, rtol=RTOL)

    def test_green_preset(self):
        # Paper section 4.2: green dust dims red and blue with D = 0.10.
        rgb = enhanced_rgb(_da([[0.0]]), _da([[0.5]]), dimming=COLOR_DIMMING["green"])
        np.testing.assert_allclose(rgb.values[0, 0, 0], 0.05 / 1.2, rtol=RTOL)
        np.testing.assert_allclose(rgb.values[0, 0, 1], 0.5 / 1.2, rtol=RTOL)
        np.testing.assert_allclose(rgb.values[0, 0, 2], 0.05 / 1.2, rtol=RTOL)

    def test_default_dimming_is_yellow(self):
        bi, cf = _da([[0.3, 0.7]]), _da([[0.2, 0.9]])
        default = enhanced_rgb(bi, cf)
        yellow = enhanced_rgb(bi, cf, dimming=COLOR_DIMMING["yellow"])
        np.testing.assert_array_equal(default.values, yellow.values)

    def test_nan_propagates(self):
        rgb = enhanced_rgb(_da([[0.5, 0.5]]), _da([[0.5, np.nan]]))
        assert np.isfinite(rgb.values[0, 0]).all()
        assert np.isnan(rgb.values[0, 1]).all()

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            enhanced_rgb(_da(np.zeros((2, 2))), _da(np.zeros((3, 3))))


class TestToUint8:
    def test_values_rounding_and_nan(self):
        rgb = xr.DataArray(
            np.array([[[0.0, 0.2, 1.0], [np.nan, 0.6, 0.75]]]),
            dims=("y", "x", "gun"),
        )
        out = to_uint8(rgb)
        assert out.dtype == np.uint8
        assert out.shape == (1, 2, 3)
        np.testing.assert_array_equal(out, [[[0, 51, 255], [0, 153, 191]]])


# --- orchestration -----------------------------------------------------------


class TestDebraImagery:
    def _inputs(self):
        area, when = object(), object()
        scene = xr.Dataset(
            {
                "refl_vis_064": _da([[10.0, 40.0], [70.0, 100.0]]),
                "bt_tir_104": _da([[250.0, 270.0], [290.0, 300.0]]),
            },
            attrs={"area": area, "start_time": when},
        )
        debra = xr.Dataset(
            {
                "cf_comb": _da([[0.0, 0.9], [0.4, 0.1]]),
                "zenith_deg": _da([[30.0, 40.0], [84.0, 120.0]]),
            }
        )
        return scene, debra, area, when

    def test_matches_composed_calls(self):
        scene, debra, area, when = self._inputs()
        out = debra_imagery(scene, debra)
        assert {"vis_bg", "ir_bg", "b_bg", "bi", "rgb"} <= set(out.data_vars)
        base = baseline_image(scene, debra["zenith_deg"])
        rgb = enhanced_rgb(base["bi"], debra["cf_comb"])
        np.testing.assert_allclose(out["bi"].values, base["bi"].values, rtol=RTOL)
        np.testing.assert_allclose(out["rgb"].values, rgb.values, rtol=RTOL)
        assert out.attrs["area"] is area
        assert out.attrs["start_time"] is when

    def test_missing_cf_raises(self):
        scene, debra, _, _ = self._inputs()
        with pytest.raises(ValueError):
            debra_imagery(scene, debra.drop_vars("cf_comb"))
        with pytest.raises(ValueError):
            debra_imagery(scene, debra.drop_vars("zenith_deg"))
