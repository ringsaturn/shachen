"""Acceptance tests for shachen.composite.composite_background.

Synthetic same-grid scene stacks with analytic expectations: per pixel the
warmest-BT10.4 *candidate* day is selected (a day is a candidate only where
all three TIR bands are finite; first wins on exact ties) and all three
background BTs are taken from that day (spectral coherence).
"""

import numpy as np
import pytest
import xarray as xr

from shachen.composite import COMPOSITE_BANDS, composite_background
from shachen.constants import Band

NAN = float("nan")


def _scene(bt86, bt104, bt123) -> xr.Dataset:
    return xr.Dataset(
        {
            "bt_tir_86": (("y", "x"), np.atleast_2d(np.asarray(bt86, dtype=float))),
            "bt_tir_104": (("y", "x"), np.atleast_2d(np.asarray(bt104, dtype=float))),
            "bt_tir_123": (("y", "x"), np.atleast_2d(np.asarray(bt123, dtype=float))),
        }
    )


def test_composite_bands_are_the_three_tir_windows():
    assert COMPOSITE_BANDS == (Band.TIR_86, Band.TIR_104, Band.TIR_123)


def test_warmest_day_selected_with_band_coherence():
    # day k has bt123 = bt104 - (1 + k), bt86 = bt104 - (2 + 2k): rsw_bg and
    # btd_bg reveal which day each pixel's *whole* band triple came from.
    day0_104 = [[290.0, 300.0], [280.0, 260.0]]
    day1_104 = [[295.0, 285.0], [270.0, 265.0]]
    scenes = [
        _scene(np.asarray(day0_104) - 2.0, day0_104, np.asarray(day0_104) - 1.0),
        _scene(np.asarray(day1_104) - 4.0, day1_104, np.asarray(day1_104) - 2.0),
    ]
    out = composite_background(scenes)

    np.testing.assert_allclose(out["bt_bg_tir_104"].values, [[295.0, 300.0], [280.0, 265.0]])
    np.testing.assert_allclose(out["rsw_bg"].values, [[-2.0, -1.0], [-1.0, -2.0]])
    np.testing.assert_allclose(out["btd_bg"].values, [[-4.0, -2.0], [-2.0, -4.0]])
    assert (out["n_valid"].values == 2).all()


def test_exact_tie_first_scene_wins():
    scenes = [
        _scene([[288.0]], [[290.0]], [[289.0]]),
        _scene([[286.0]], [[290.0]], [[287.0]]),
    ]
    out = composite_background(scenes)
    np.testing.assert_allclose(out["rsw_bg"].values, [[-1.0]])
    np.testing.assert_allclose(out["btd_bg"].values, [[-2.0]])


def test_day_with_any_nan_band_is_not_a_candidate():
    # day1 is warmer in bt104 but has NaN bt86 there -> day0 must be selected.
    scenes = [
        _scene([[288.0]], [[290.0]], [[289.0]]),
        _scene([[NAN]], [[295.0]], [[294.0]]),
    ]
    out = composite_background(scenes)
    np.testing.assert_allclose(out["bt_bg_tir_104"].values, [[290.0]])
    np.testing.assert_allclose(out["bt_bg_tir_86"].values, [[288.0]])
    np.testing.assert_allclose(out["bt_bg_tir_123"].values, [[289.0]])
    assert (out["n_valid"].values == 1).all()


def test_zero_candidate_pixels_are_nan():
    # pixel 0: candidate in day0 only (day1 bt86 is NaN there); pixel 1: NaN
    # in a band of every day.
    scenes = [
        _scene([[288.0, 279.0]], [[290.0, NAN]], [[289.0, 280.0]]),
        _scene([[NAN, 288.0]], [[285.0, 290.0]], [[284.0, NAN]]),
    ]
    out = composite_background(scenes)
    np.testing.assert_allclose(out["bt_bg_tir_104"].values[0, 0], 290.0)
    assert out["n_valid"].values[0, 0] == 1
    for name in ("bt_bg_tir_86", "bt_bg_tir_104", "bt_bg_tir_123", "rsw_bg", "btd_bg"):
        assert np.isnan(out[name].values[0, 1]), name
    assert out["n_valid"].values[0, 1] == 0


def test_single_scene_is_identity():
    scene = _scene([[286.0, 284.0]], [[290.0, 288.0]], [[288.5, 286.5]])
    out = composite_background([scene])
    np.testing.assert_allclose(out["bt_bg_tir_104"].values, scene["bt_tir_104"].values)
    np.testing.assert_allclose(out["rsw_bg"].values, [[-1.5, -1.5]])
    np.testing.assert_allclose(out["btd_bg"].values, [[-4.0, -4.0]])
    assert (out["n_valid"].values == 1).all()


def test_output_contract_names_attrs_dtype():
    scenes = [
        _scene([[288.0]], [[290.0]], [[289.0]]),
        _scene([[287.0]], [[289.0]], [[288.0]]),
    ]
    out = composite_background(scenes)
    assert {
        "bt_bg_tir_86",
        "bt_bg_tir_104",
        "bt_bg_tir_123",
        "rsw_bg",
        "btd_bg",
        "n_valid",
    } <= set(out.data_vars)
    for name in ("bt_bg_tir_86", "bt_bg_tir_104", "bt_bg_tir_123"):
        assert out[name].attrs["units"] == "K"
    assert np.issubdtype(out["n_valid"].dtype, np.integer)
    assert out.attrs["n_scenes"] == 2


def test_empty_sequence_raises():
    with pytest.raises(ValueError):
        composite_background([])


def test_missing_variable_raises():
    good = _scene([[288.0]], [[290.0]], [[289.0]])
    bad = good.drop_vars("bt_tir_86")
    with pytest.raises(ValueError):
        composite_background([good, bad])


def test_shape_mismatch_raises():
    a = _scene([[288.0]], [[290.0]], [[289.0]])
    b = _scene([[288.0, 287.0]], [[290.0, 289.0]], [[289.0, 288.0]])
    with pytest.raises(ValueError):
        composite_background([a, b])
