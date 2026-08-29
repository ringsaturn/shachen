"""Acceptance tests for shachen.dustrgb (classic EUMETSAT Dust RGB baseline).

All fields are synthetic and every expectation is an exact algebraic value of
the frozen recipe (constants.DUST_RGB: R = BT12.3-BT10.4 on [-4, 2] K,
G = BT11.2-BT8.6 on [0, 15] K gamma 2.5, B = BT10.4 on [261, 289] K), so
tolerances are tight.
"""

import numpy as np
import pytest
import xarray as xr

from shachen.constants import DUST_RGB, DUST_RGB_ABI
from shachen.dustrgb import dust_rgb
from shachen.pipeline import run_dust_rgb

RTOL = 1e-9


def _da(values) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=float), dims=("y", "x"))


def _scene(bt86, bt104, bt112, bt123) -> xr.Dataset:
    return xr.Dataset(
        {
            "bt_tir_86": _da(bt86),
            "bt_tir_104": _da(bt104),
            "bt_tir_112": _da(bt112),
            "bt_tir_123": _da(bt123),
        }
    )


class TestRecipe:
    def test_red_split_window_stretch(self):
        # BT12.3 - BT10.4 = -4 / -1 / +2 K -> 0 / 0.5 / 1 (gamma 1).
        bt104 = [[280.0, 280.0, 280.0]]
        scene = _scene(bt104, bt104, bt104, [[276.0, 279.0, 282.0]])
        out = dust_rgb(scene)
        np.testing.assert_allclose(out.sel(gun="r").values, [[0.0, 0.5, 1.0]], rtol=RTOL)

    def test_green_uses_112_with_gamma(self):
        # BT11.2 - BT8.6 = 0 / 7.5 / 15 K -> 0 / 0.5**(1/2.5) / 1.
        bt104 = [[280.0, 280.0, 280.0]]
        scene = _scene([[280.0, 280.0, 280.0]], bt104, [[280.0, 287.5, 295.0]], bt104)
        out = dust_rgb(scene)
        expected = [[0.0, 0.5 ** (1.0 / DUST_RGB.green_gamma), 1.0]]
        np.testing.assert_allclose(out.sel(gun="g").values, expected, rtol=RTOL)

    def test_blue_bt104_stretch(self):
        # BT10.4 = 261 / 275 / 289 K -> 0 / 0.5 / 1.
        bt104 = [[261.0, 275.0, 289.0]]
        scene = _scene(bt104, bt104, bt104, bt104)
        out = dust_rgb(scene)
        np.testing.assert_allclose(out.sel(gun="b").values, [[0.0, 0.5, 1.0]], rtol=RTOL)

    def test_out_of_range_clips(self):
        # Way outside every stretch: guns saturate at 0/1, never leave [0, 1].
        scene = _scene([[320.0]], [[200.0]], [[200.0]], [[320.0]])
        out = dust_rgb(scene)
        np.testing.assert_allclose(out.values, [[[1.0, 0.0, 0.0]]], rtol=RTOL)


class TestLayout:
    def test_dims_and_gun_coord_match_enhanced_rgb(self):
        bt = [[270.0, 280.0], [275.0, 285.0]]
        out = dust_rgb(_scene(bt, bt, bt, bt))
        assert out.dims == ("y", "x", "gun")
        assert list(out.coords["gun"].values) == ["r", "g", "b"]

    def test_nan_propagates(self):
        bt = [[280.0, np.nan]]
        out = dust_rgb(_scene(bt, bt, bt, bt))
        assert np.isfinite(out.values[0, 0]).all()
        assert np.isnan(out.values[0, 1]).all()

    def test_missing_variable_raises(self):
        scene = _scene([[280.0]], [[280.0]], [[280.0]], [[280.0]]).drop_vars("bt_tir_112")
        with pytest.raises(ValueError, match="bt_tir_112"):
            dust_rgb(scene)


class TestPipelineEntry:
    """run_dust_rgb is the pipeline-level entry, mirroring run_debra."""

    def test_wraps_the_composite_with_scene_attrs(self):
        bt = [[270.0, 280.0], [275.0, 285.0]]
        scene = _scene(bt, bt, bt, bt)
        scene.attrs.update(area="fake-area", start_time="2017-03-23T21:00")

        out = run_dust_rgb(scene)
        assert list(out.data_vars) == ["dust_rgb"]
        assert out.attrs["area"] == "fake-area"
        assert out.attrs["start_time"] == "2017-03-23T21:00"

    def test_debra_only_scene_raises(self):
        # A scene loaded with roles=DEBRA_BANDS has no 11.2 um band.
        bt = [[280.0]]
        scene = _scene(bt, bt, bt, bt).drop_vars("bt_tir_112")
        with pytest.raises(ValueError, match="bt_tir_112"):
            run_dust_rgb(scene)


class TestPerSensorStretches:
    """The scheme is retuned per imager; the reader attr picks the set."""

    @staticmethod
    def _typical_scene(reader=None):
        # Clear-sky-ish land: the two sets disagree most on the red gun.
        scene = _scene([[275.0]], [[288.0]], [[286.0]], [[285.0]])
        if reader is not None:
            scene.attrs["reader"] = reader
        return scene

    @pytest.mark.parametrize(
        ("reader", "expected"),
        [("abi_l1b", DUST_RGB_ABI), ("ahi_hsd", DUST_RGB)],
    )
    def test_reader_selects_the_operational_set(self, reader, expected):
        out = run_dust_rgb(self._typical_scene(reader))
        np.testing.assert_allclose(
            out["dust_rgb"].values,
            dust_rgb(self._typical_scene(), expected).values,
            rtol=RTOL,
        )

    @pytest.mark.parametrize("reader", [None, "seviri_l1b5"])
    def test_unknown_or_absent_reader_falls_back_to_seviri(self, reader):
        out = run_dust_rgb(self._typical_scene(reader))
        np.testing.assert_allclose(
            out["dust_rgb"].values,
            dust_rgb(self._typical_scene(), DUST_RGB).values,
            rtol=RTOL,
        )

    def test_explicit_constants_override_the_reader(self):
        out = run_dust_rgb(self._typical_scene("abi_l1b"), DUST_RGB)
        np.testing.assert_allclose(
            out["dust_rgb"].values,
            dust_rgb(self._typical_scene(), DUST_RGB).values,
            rtol=RTOL,
        )

    def test_the_two_sets_actually_differ(self):
        # Guards the whole point: if these ever coincide the dispatch is dead
        # code and the per-sensor rule has silently stopped meaning anything.
        scene = self._typical_scene()
        assert not np.allclose(
            dust_rgb(scene, DUST_RGB).values, dust_rgb(scene, DUST_RGB_ABI).values
        )
