"""MERRA-2 cache loading: time interpolation between the hourly means."""

import datetime as dt

import numpy as np
import pytest
import xarray as xr

from shachen.io.merra import load_skin_temperature, load_surface_meteorology


@pytest.fixture
def day_file(tmp_path):
    """A one-day TS + met cache stamped like MERRA-2: 00:30, 01:30, ... 23:30.

    TS rises 1 K per hour from 280 K so an interpolated value names its own
    time, and the met variables carry the matchup's covariate set.
    """
    times = [np.datetime64("2021-03-15T00:30", "ns") + np.timedelta64(h, "h") for h in range(24)]
    lat = np.array([30.0, 40.0])
    lon = np.array([100.0, 110.0])
    shape = (len(times), lat.size, lon.size)
    hours = np.arange(len(times), dtype=float).reshape(-1, 1, 1)
    dataset = xr.Dataset(
        {
            "TS": (("time", "lat", "lon"), np.broadcast_to(280.0 + hours, shape).copy()),
            "T2M": (("time", "lat", "lon"), np.broadcast_to(270.0 + hours, shape).copy()),
            "PBLH": (("time", "lat", "lon"), np.full(shape, 500.0)),
        },
        coords={"time": times, "lat": lat, "lon": lon},
    )
    path = tmp_path / "merra2_ts_20210315.nc"
    dataset.to_netcdf(path)
    return path


def test_skin_temperature_interpolates_between_stamps(day_file):
    """02:00 sits halfway between the 01:30 and 02:30 hourly means."""
    ts = load_skin_temperature(day_file, dt.datetime(2021, 3, 15, 2, 0))
    assert ts.values == pytest.approx(281.5)
    assert ts.attrs["units"] == "K"


def test_skin_temperature_clamps_at_the_start_of_the_day(day_file):
    """00:00 UTC is before the first stamp: hold it, never extrapolate to NaN.

    Regression: plain interp returned all-NaN there, which silently emptied
    the background, the cloud mask and every dust test for 00 UTC scans.
    """
    ts = load_skin_temperature(day_file, dt.datetime(2021, 3, 15, 0, 0))
    assert np.isfinite(ts.values).all()
    assert ts.values == pytest.approx(280.0)  # the 00:30 mean


def test_skin_temperature_clamps_at_the_end_of_the_day(day_file):
    """A scan after the last stamp (23:30) holds that hourly mean too."""
    ts = load_skin_temperature(day_file, dt.datetime(2021, 3, 15, 23, 50))
    assert ts.values == pytest.approx(303.0)


def test_surface_meteorology_clamps_the_same_way(day_file):
    """The PM10 covariates share the fix; NaN there empties the retrieval."""
    met = load_surface_meteorology(day_file, dt.datetime(2021, 3, 15, 0, 0))
    assert met["T2M"].values == pytest.approx(270.0)
    assert met["PBLH"].values == pytest.approx(500.0)
    assert np.isfinite(met["T2M"].values).all()
