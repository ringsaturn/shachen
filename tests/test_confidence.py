"""Acceptance tests for debra.confidence (Eqs. 16-22, erratum forms of 21-22).

Pure float64 arithmetic -> rtol=1e-6, except assertions involving the
trigonometric blend weights, which use precomputed goldens at rtol=1e-4:
b_ngt_trm(97.5 deg) = ((cos 97.5 - cos 105) / (cos 90 - cos 105))^1.5
              = 0.348987 (NOAA-independent hand computation).
Inputs derive from DEFAULTS.confidence where bounds are involved.
"""

import dataclasses

import numpy as np
import pytest
import xarray as xr

from shachen.confidence import confidence, confidence_norm, confidence_raw
from shachen.constants import DEFAULTS, Bounds, ConfidenceConstants

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


def _norm(raw, bounds=None):
    """Eq. 19 with the day interval by default."""
    bounds = bounds or C.cf_norm_day
    return np.clip((raw - bounds.min) / (bounds.max - bounds.min), 0.0, 1.0)


def _norm_ngt(raw):
    return _norm(raw, C.cf_norm_ngt)


def _norm_trm(raw):
    """Eq. 19 on the terminator's own interval (Eq. 17 ceiling 2.5)."""
    return _norm(raw, C.cf_norm_trm)


def _norm_interp(raw, weight):
    """Eq. 19 with the interval interpolated -- the 0.3.0 terminator."""
    lower = weight * C.cf_norm_day.min + (1.0 - weight) * C.cf_norm_ngt.min
    upper = weight * C.cf_norm_day.max + (1.0 - weight) * C.cf_norm_ngt.max
    return np.clip((raw - lower) / (upper - lower), 0.0, 1.0)


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
    # DT sum at the midpoint of the day interval -> cf_day = 0.5 exactly.
    mid = (C.cf_norm_day.min + C.cf_norm_day.max) / 2.0
    out = confidence(*_inputs(mid - 1.0, 0.5, 0.5, cm_day=0.0, zen=30.0))
    np.testing.assert_allclose(out["cf_comb"].values, 0.5, rtol=1e-6)


def test_night_uses_max_and_night_mask():
    # Eq. 18: max(DT1, DT2) + 0.5*DT3, suppressed by cm_norm_ngt (not _day).
    out = confidence(*_inputs(0.2, 0.8, 0.4, cm_day=1.0, cm_ngt=0.0, zen=120.0))
    raw = max(0.2, 0.8) + C.dt3_weight_ngt * 0.4  # = 1.0
    np.testing.assert_allclose(out["b_trm_day"].values, 0.0, atol=1e-12)
    np.testing.assert_allclose(out["b_ngt_trm"].values, 0.0, atol=1e-12)
    np.testing.assert_allclose(out["cf_comb"].values, _norm_ngt(raw), rtol=1e-6)


def test_terminator_selects_cf_trm():
    # At zenith 90: b_trm_day = 0, b_ngt_trm = 1 -> cf_comb == cf_trm (Eq. 17),
    # normalized on the terminator's own interval.
    out = confidence(*_inputs(0.6, 0.2, 0.8, cm_day=0.0, zen=90.0))
    raw = 0.6 + 0.2 + C.dt3_weight_trm * 0.8  # = 1.2
    np.testing.assert_allclose(out["cf_trm"].values, _norm_trm(raw), rtol=1e-6)
    np.testing.assert_allclose(out["cf_comb"].values, out["cf_trm"].values, rtol=1e-6)


def test_blend_weight_golden_and_composition():
    # Eq. 22 (erratum) at zenith 97.5 deg with distinct cf_day/trm/ngt.
    out = confidence(*_inputs(0.6, 0.2, 0.8, cm_day=0.0, cm_ngt=0.5, zen=97.5))
    np.testing.assert_allclose(out["b_ngt_trm"].values, B_NGT_TRM_97_5, rtol=1e-4)
    np.testing.assert_allclose(out["b_trm_day"].values, 0.0, atol=1e-12)
    cf_trm = _norm_trm(0.6 + 0.2 + C.dt3_weight_trm * 0.8)
    cf_ngt = _norm_ngt((max(0.6, 0.2) + C.dt3_weight_ngt * 0.8) * (1.0 - 0.5))
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


def test_night_branch_reaches_one():
    # Plan 008 step 2: with its own Eq. 19 interval the night branch spans the
    # whole of [0, 1]. Eq. 18 at its ceiling -- max(DT1, DT2) = 1 and DT3 = 1,
    # i.e. raw 1.5 -- must saturate, not stall near 0.5 as the shared interval
    # made it do ((1.5 - 0.25) / (2.50 - 0.25) = 0.556).
    out = confidence(*_inputs(1.0, 1.0, 1.0, cm_ngt=0.0, zen=140.0))
    np.testing.assert_allclose(out["cf_ngt"].values, 1.0, rtol=1e-6)
    np.testing.assert_allclose(out["cf_comb"].values, 1.0, rtol=1e-6)


def test_every_branch_agrees_on_a_shared_signal():
    # The point of the split: the same normalized dust signal reads the same
    # in every branch. DT1 = DT2 = DT3 = v gives raws of 3v, 2.5v and 1.5v,
    # and each interval is scaled to its own ceiling, so all three normalize
    # to the same number -- including at 90 deg, where CF_comb is CF_trm.
    for value in (0.2, 0.5, 0.9):
        day = confidence(*_inputs(value, value, value, zen=20.0))
        terminator = confidence(*_inputs(value, value, value, zen=90.0))
        night = confidence(*_inputs(value, value, value, zen=140.0))
        np.testing.assert_allclose(night["cf_comb"].values, day["cf_comb"].values, rtol=1e-6)
        np.testing.assert_allclose(terminator["cf_comb"].values, day["cf_comb"].values, rtol=1e-6)


def test_terminator_interval_can_be_switched_back_to_interpolation():
    # cf_norm_trm=None is 0.3.0: no interval of its own, the Eq. 20 weight
    # interpolating the day and night ones. Kept so the change the 42-day
    # evaluation scored is reproducible from the same code.
    legacy = dataclasses.replace(C, cf_norm_trm=None)
    out = confidence(*_inputs(0.6, 0.2, 0.8, cm_day=0.0, zen=97.5), constants=legacy)
    raw = 0.6 + 0.2 + C.dt3_weight_trm * 0.8
    np.testing.assert_allclose(out["cf_trm"].values, _norm_interp(raw, B_NGT_TRM_97_5), rtol=1e-4)


def test_the_terminator_read_low_before_its_own_interval():
    # What the 0.3.0 terminator cost: on the day side of 90 deg CF_trm was
    # normalized on the day interval, whose ceiling is 3.0, while Eq. 17 only
    # reaches 2.5 -- so the same dust read lower at dusk than it did at noon.
    legacy = dataclasses.replace(C, cf_norm_trm=None)
    args = _inputs(0.5, 0.5, 0.5, zen=90.0)
    assert confidence(*args, constants=legacy)["cf_trm"].values.max() < (
        confidence(*args)["cf_trm"].values.min()
    )


def test_cf_norm_alias_restores_the_single_interval():
    # Backward compatibility: passing cf_norm sets both intervals, which is
    # the pre-split behaviour -- the night ceiling drops back to 0.556.
    legacy = ConfidenceConstants(cf_norm=Bounds(0.25, 2.50))
    assert legacy.cf_norm_day == legacy.cf_norm_ngt == Bounds(0.25, 2.50)
    out = confidence(*_inputs(1.0, 1.0, 1.0, zen=140.0), constants=legacy)
    np.testing.assert_allclose(out["cf_ngt"].values, (1.5 - 0.25) / 2.25, rtol=1e-6)


def test_cf_norm_does_not_read_back_as_an_interval():
    # Not a silent alias: it is an InitVar, so code that used to read
    # constants.cf_norm.min raises instead of normalizing with wrong bounds.
    assert DEFAULTS.confidence.cf_norm is None
    with pytest.raises(AttributeError):
        _ = DEFAULTS.confidence.cf_norm.min


def test_replace_still_works():
    # The InitVar stays declared on the class so dataclasses.replace, which
    # reads every init field off the instance, does not trip over it.
    tuned = dataclasses.replace(DEFAULTS.confidence, cf_norm_ngt=Bounds(0.1, 1.0))
    assert tuned.cf_norm_ngt == Bounds(0.1, 1.0)
    assert tuned.cf_norm_day == DEFAULTS.confidence.cf_norm_day


def test_confidence_raw_matches_the_normalized_branches():
    # Plan 008 step 4: the raw Eqs. 16-18 sums are what an offline interval
    # search re-normalizes, so they have to be the very sums confidence()
    # normalizes -- not a second copy of the equations.
    tests, cloud, zenith = _inputs(0.6, 0.2, 0.8, cm_day=0.3, cm_ngt=0.5, zen=97.5)
    raw = confidence_raw(tests, cloud)
    out = confidence(tests, cloud, zenith)
    np.testing.assert_allclose(raw["cf_day_raw"].values, (0.6 + 0.2 + 0.8) * (1.0 - 0.3), rtol=1e-6)
    np.testing.assert_allclose(
        raw["cf_ngt_raw"].values, (0.6 + C.dt3_weight_ngt * 0.8) * (1.0 - 0.5), rtol=1e-6
    )
    np.testing.assert_allclose(out["cf_day"].values, _norm(raw["cf_day_raw"].values), rtol=1e-6)
    np.testing.assert_allclose(out["cf_ngt"].values, _norm_ngt(raw["cf_ngt_raw"].values), rtol=1e-6)
    np.testing.assert_allclose(out["cf_trm"].values, _norm_trm(raw["cf_trm_raw"].values), rtol=1e-6)


def test_confidence_raw_does_not_depend_on_the_intervals():
    # Eq. 19 is downstream of these sums: re-tuning the intervals must not
    # move them, which is what makes one sampled archive enough for a grid.
    tests, cloud, _ = _inputs(0.6, 0.2, 0.8, cm_day=0.3, cm_ngt=0.5)
    other = dataclasses.replace(C, cf_norm_day=Bounds(0.0, 1.0), cf_norm_ngt=Bounds(0.05, 0.9))
    for name, values in confidence_raw(tests, cloud, other).items():
        np.testing.assert_allclose(values.values, confidence_raw(tests, cloud)[name].values)


def test_confidence_raw_requires_its_inputs():
    tests, cloud, _ = _inputs(0.5, 0.5, 0.5)
    with pytest.raises(ValueError):
        confidence_raw(tests.drop_vars("dt3"), cloud)
    with pytest.raises(ValueError):
        confidence_raw(tests, cloud.drop_vars("cm_norm_ngt"))


def test_confidence_norm_reproduces_confidence():
    # The two halves compose back into the whole: Eqs. 16-18 then Eqs. 19-22.
    for zen in (30.0, 90.0, 97.5, 140.0):
        tests, cloud, zenith = _inputs(0.6, 0.2, 0.8, cm_day=0.3, cm_ngt=0.5, zen=zen)
        whole = confidence(tests, cloud, zenith)
        halves = confidence_norm(confidence_raw(tests, cloud), zenith)
        for name in EXPECTED_VARS:
            np.testing.assert_allclose(halves[name].values, whole[name].values, rtol=1e-12)


def test_confidence_norm_retunes_without_the_dust_tests():
    # What plan 008 step 4 needs: score a different Eq. 19 night interval on
    # stored raw sums, with no access to DT1-DT3 or the L1b behind them.
    tests, cloud, zenith = _inputs(0.6, 0.2, 0.8, cm_ngt=0.0, zen=140.0)
    raw = confidence_raw(tests, cloud)
    tuned = dataclasses.replace(C, cf_norm_ngt=Bounds(0.05, 1.0))
    out = confidence_norm(raw, zenith, tuned)
    expected = np.clip((raw["cf_ngt_raw"].values - 0.05) / (1.0 - 0.05), 0.0, 1.0)
    np.testing.assert_allclose(out["cf_comb"].values, expected, rtol=1e-6)


def test_confidence_norm_requires_the_raw_sums():
    tests, cloud, zenith = _inputs(0.5, 0.5, 0.5)
    with pytest.raises(ValueError):
        confidence_norm(confidence_raw(tests, cloud).drop_vars("cf_trm_raw"), zenith)
