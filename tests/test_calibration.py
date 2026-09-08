"""CF_cal: the monotone tables and the field they produce.

No ground truth here -- fitting the knots needs stations and lives outside
this package. What is testable is that the table is a monotone map with the
end conditions it claims, that CF_cal is CF_day wherever Eq. 22 reads the day
branch, and that the shipped table keeps the property a re-tuned Eq. 19
interval could not: a saturated night branch still reads 1.0.
"""

import numpy as np
import pytest
import xarray as xr

from shachen.calibration import EAST_ASIA_AHI, Calibration, CalibrationTable, calibrate
from shachen.confidence import confidence, confidence_raw
from shachen.constants import DEFAULTS

C = DEFAULTS.confidence
SHAPE = (2, 2)

TABLE = CalibrationTable(knots=((0.2, 0.0), (0.6, 0.5), (1.5, 1.0)), calibrated_to=0.5)


def _raw(day=1.0, trm=1.0, ngt=1.0):
    return xr.Dataset(
        {
            "cf_day_raw": xr.DataArray(np.full(SHAPE, day), dims=("y", "x")),
            "cf_trm_raw": xr.DataArray(np.full(SHAPE, trm), dims=("y", "x")),
            "cf_ngt_raw": xr.DataArray(np.full(SHAPE, ngt), dims=("y", "x")),
        }
    )


class TestTable:
    def test_interpolates_between_knots(self):
        np.testing.assert_allclose(TABLE(0.4), 0.25, rtol=1e-12)
        np.testing.assert_allclose(TABLE(1.05), 0.75, rtol=1e-12)

    def test_clamps_outside_the_knots(self):
        assert TABLE(0.0) == 0.0
        assert TABLE(-3.0) == 0.0
        assert TABLE(99.0) == 1.0

    def test_nan_propagates(self):
        assert np.isnan(TABLE(np.array([np.nan, 0.6]))[0])
        assert TABLE(np.array([np.nan, 0.6]))[1] == pytest.approx(0.5)

    def test_dataarray_keeps_its_coords(self):
        raw = xr.DataArray(np.full(SHAPE, 0.6), dims=("y", "x"), coords={"y": [0, 1], "x": [0, 1]})
        out = TABLE(raw)
        assert isinstance(out, xr.DataArray)
        assert out.dims == ("y", "x")
        np.testing.assert_allclose(out.values, 0.5)

    @pytest.mark.parametrize(
        "knots",
        [
            (((0.2, 0.0),)),  # one knot is not a map
            (((0.6, 0.0), (0.2, 1.0))),  # raw must increase
            (((0.2, 0.0), (0.2, 1.0))),  # strictly: no vertical jumps
            (((0.2, 0.0), (0.6, 0.7), (1.5, 0.5))),  # cf must not decrease
            (((0.2, 0.1), (1.5, 1.0))),  # must start at 0
            (((0.2, 0.0), (1.5, 0.9))),  # must end at 1
        ],
    )
    def test_rejects_a_table_that_is_not_a_monotone_map_onto_0_1(self, knots):
        with pytest.raises(ValueError):
            CalibrationTable(knots=knots)


class TestCalibrate:
    def test_day_is_the_reference_and_is_left_alone(self):
        raw, zenith = _raw(day=1.45), xr.DataArray(np.full(SHAPE, 20.0), dims=("y", "x"))
        out = calibrate(raw, zenith, Calibration(night=TABLE, terminator=TABLE), C)
        np.testing.assert_allclose(out["b_trm_day"].values, 1.0, rtol=1e-9)
        expected = (1.45 - C.cf_norm_day.min) / (C.cf_norm_day.max - C.cf_norm_day.min)
        np.testing.assert_allclose(out["cf_cal"].values, expected, rtol=1e-9)
        np.testing.assert_allclose(out["cf_cal"].values, out["cf_day"].values, rtol=1e-12)

    def test_night_reads_the_table(self):
        raw, zenith = _raw(ngt=0.6), xr.DataArray(np.full(SHAPE, 140.0), dims=("y", "x"))
        out = calibrate(raw, zenith, Calibration(night=TABLE, terminator=TABLE), C)
        np.testing.assert_allclose(out["cf_cal"].values, 0.5, rtol=1e-12)

    def test_blend_is_eq_22_on_the_calibrated_branches(self):
        raw = _raw(day=1.0, trm=0.6, ngt=1.5)
        zenith = xr.DataArray(np.full(SHAPE, 97.5), dims=("y", "x"))
        out = calibrate(raw, zenith, Calibration(night=TABLE, terminator=TABLE), C)
        weight = out["b_ngt_trm"].values
        expected = weight * 0.5 + (1.0 - weight) * 1.0
        np.testing.assert_allclose(out["cf_cal"].values, expected, rtol=1e-12)

    def test_missing_raw_variable_raises(self):
        raw = _raw().drop_vars("cf_ngt_raw")
        with pytest.raises(ValueError):
            calibrate(raw, xr.DataArray(np.full(SHAPE, 30.0), dims=("y", "x")))

    def test_cf_comb_is_not_touched(self):
        # The product decision is that both fields exist: CF_comb keeps the
        # physical-fraction meaning, CF_cal carries the precision one.
        tests = xr.Dataset(
            {
                name: xr.DataArray(np.full(SHAPE, 0.7), dims=("y", "x"))
                for name in ("dt1", "dt2", "dt3")
            }
        )
        cloud = xr.Dataset(
            {
                name: xr.DataArray(np.zeros(SHAPE), dims=("y", "x"))
                for name in ("cm_norm_day", "cm_norm_ngt")
            }
        )
        zenith = xr.DataArray(np.full(SHAPE, 140.0), dims=("y", "x"))
        before = confidence(tests, cloud, zenith).copy(deep=True)
        calibrate(confidence_raw(tests, cloud), zenith)
        after = confidence(tests, cloud, zenith)
        for name in before.data_vars:
            np.testing.assert_array_equal(before[name].values, after[name].values)


class TestShippedTable:
    def test_a_saturated_night_branch_still_reads_one(self):
        # What the linear alternative could not do: matching day precision by
        # moving the Eq. 19 interval capped CF_ngt near 0.55. The table pins
        # the working points *and* keeps Eq. 18's ceiling of 1.5 at 1.0.
        assert EAST_ASIA_AHI.night(1.5) == pytest.approx(1.0)
        assert EAST_ASIA_AHI.terminator(2.5) == pytest.approx(1.0)

    def test_the_calibrated_range_is_declared(self):
        # Above 0.45 night precision has nothing to match, so the table says
        # so rather than implying a precision it cannot deliver.
        assert EAST_ASIA_AHI.night.calibrated_to == 0.45
        assert EAST_ASIA_AHI.night(0.9238) == pytest.approx(0.45)

    def test_it_is_monotone_over_the_whole_raw_range(self):
        raw = np.linspace(0.0, 1.6, 400)
        assert np.all(np.diff(EAST_ASIA_AHI.night(raw)) >= 0.0)
