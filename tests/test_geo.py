"""Acceptance tests for shachen.geo: ancillary regridding and the land mask.

regrid_latlon is the bilinear bridge from the MERRA-2 / CAMEL lat/lon grids to
the satellite grid; land_mask feeds DT3's land/ocean surface shift (Eq. 15).
Tolerance: bilinear interpolation reproduces an affine field f(lat, lon)
exactly, so atol=1e-6 (pure float64 arithmetic headroom).
"""

import numpy as np
import xarray as xr

from shachen.geo import land_mask, regrid_latlon


def _affine_field(lats, lons, lat_name="lat", lon_name="lon"):
    field = 2.0 * lats[:, None] + 0.5 * lons[None, :]
    return xr.DataArray(field, coords={lat_name: lats, lon_name: lons}, dims=(lat_name, lon_name))


def test_regrid_affine_field_exact(latlon_area):
    src = _affine_field(np.arange(25.0, 41.0), np.arange(-110.0, -94.0))
    out = regrid_latlon(src, latlon_area)
    lons, lats = latlon_area.get_lonlats()
    assert out.dims == ("y", "x")
    assert out.shape == (latlon_area.height, latlon_area.width)
    np.testing.assert_allclose(out.values, 2.0 * lats + 0.5 * lons, atol=1e-6)


def test_regrid_accepts_latitude_longitude_dim_names(latlon_area):
    src = _affine_field(
        np.arange(25.0, 41.0),
        np.arange(-110.0, -94.0),
        lat_name="latitude",
        lon_name="longitude",
    )
    out = regrid_latlon(src, latlon_area)
    lons, lats = latlon_area.get_lonlats()
    np.testing.assert_allclose(out.values, 2.0 * lats + 0.5 * lons, atol=1e-6)


def test_regrid_dataset_regrids_every_variable(latlon_area):
    lats, lons = np.arange(25.0, 41.0), np.arange(-110.0, -94.0)
    ds = xr.Dataset(
        {
            "a": _affine_field(lats, lons),
            "b": 3.0 * _affine_field(lats, lons),
        }
    )
    out = regrid_latlon(ds, latlon_area)
    tlons, tlats = latlon_area.get_lonlats()
    np.testing.assert_allclose(out["a"].values, 2.0 * tlats + 0.5 * tlons, atol=1e-6)
    np.testing.assert_allclose(out["b"].values, 3.0 * (2.0 * tlats + 0.5 * tlons), atol=1e-6)


def test_regrid_outside_source_extent_is_nan(latlon_area):
    # Source only covers 32-40N; target pixel centers south of 32N have no data.
    src = _affine_field(np.arange(32.0, 41.0), np.arange(-110.0, -94.0))
    out = regrid_latlon(src, latlon_area)
    _, lats = latlon_area.get_lonlats()
    assert np.isnan(out.values[lats < 32.0]).all()
    assert np.isfinite(out.values[lats > 32.0]).all()


def test_regrid_offdisk_is_nan(geos_area):
    src = _affine_field(np.arange(-90.0, 91.0, 5.0), np.arange(-180.0, 181.0, 5.0))
    out = regrid_latlon(src, geos_area)
    assert np.isnan(out.values[0, 0])  # corner is off the Earth disk
    assert np.isfinite(out.values[25, 25])  # disk center


def test_land_mask_all_land(latlon_area):
    mask = land_mask(latlon_area)
    assert mask.dims == ("y", "x")
    assert mask.shape == (latlon_area.height, latlon_area.width)
    assert mask.dtype == bool
    assert mask.values.all()


def test_land_mask_all_ocean(make_latlon_area):
    gulf = make_latlon_area(-92.0, 24.0, -90.0, 26.0, 8, 8)  # Gulf of Mexico
    assert not land_mask(gulf).values.any()


def test_land_mask_offdisk_is_false(geos_area):
    mask = land_mask(geos_area)
    assert not mask.values[0, 0]  # off-disk corner
    assert mask.values.any()  # the disk does contain land
