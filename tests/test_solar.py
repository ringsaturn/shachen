"""Acceptance tests for shachen.solar.solar_zenith.

Golden values were computed with the NOAA solar position algorithm (the NOAA
GML spreadsheet equations), independently of the pyorbital implementation
path. Tolerance abs=0.5 deg: the day/terminator/night blend spans 30 deg of
zenith (Eqs. 20-21), so half a degree has no algorithmic effect.
"""

import datetime as dt

import numpy as np
import pytest

from shachen.solar import solar_zenith

GOLDEN = [
    # (lat, lon, when, zenith_deg per NOAA SPA)
    (35.0, -100.0, dt.datetime(2017, 3, 23, 18, 0), 35.3823),  # local late morning
    (35.0, -100.0, dt.datetime(2017, 3, 23, 9, 0), 132.2486),  # deep night
    (0.0, 0.0, dt.datetime(2017, 3, 22, 12, 0), 1.8834),  # near-subsolar, equinox
]


@pytest.mark.parametrize("lat,lon,when,expected", GOLDEN)
def test_golden_points(make_latlon_area, lat, lon, when, expected):
    area = make_latlon_area(lon - 0.5, lat - 0.5, lon + 0.5, lat + 0.5, 1, 1)
    zen = solar_zenith(area, when)
    assert zen.dims == ("y", "x")
    assert zen.shape == (1, 1)
    assert zen.values[0, 0] == pytest.approx(expected, abs=0.5)


def test_grid_shape_and_range(latlon_area):
    zen = solar_zenith(latlon_area, dt.datetime(2017, 3, 23, 18, 0))
    assert zen.dims == ("y", "x")
    assert zen.shape == (latlon_area.height, latlon_area.width)
    assert ((zen.values >= 0.0) & (zen.values <= 180.0)).all()


def test_offdisk_is_nan(geos_area):
    zen = solar_zenith(geos_area, dt.datetime(2017, 3, 23, 18, 0))
    assert np.isnan(zen.values[0, 0])  # off-disk corner
    assert np.isfinite(zen.values[25, 25])  # disk center
