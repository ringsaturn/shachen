"""Grid helpers: regrid lat/lon ancillary fields to the satellite grid, land mask.

MERRA-2 skin temperature (0.5 x 0.625 deg) and the CAMEL emissivity climatology
(0.05 deg) live on regular lat/lon grids; every DEBRA computation happens on the
sensor's 2-km fixed grid, so ancillary data is interpolated onto the scene's
pyresample area (never the reverse). The land mask feeds DT3's land/ocean
surface shift (Eq. 15).
"""

import numpy as np
import xarray as xr
from global_land_mask import globe
from pyresample.geometry import AreaDefinition

_LAT_NAMES = ("lat", "latitude")
_LON_NAMES = ("lon", "longitude")


def _find_dim(source: xr.DataArray | xr.Dataset, candidates: tuple[str, ...]) -> str:
    for name in candidates:
        if name in source.dims:
            return name
    raise KeyError(
        f"source has none of the expected dimensions {candidates}; got {tuple(source.dims)}"
    )


def _target_lonlats(area: AreaDefinition) -> tuple[np.ndarray, np.ndarray]:
    """Target lon/lat as float64 arrays with off-disk pixels set to NaN."""
    lons, lats = area.get_lonlats()
    lons = np.ma.filled(np.asarray(lons, dtype=np.float64), np.nan)
    lats = np.ma.filled(np.asarray(lats, dtype=np.float64), np.nan)
    offdisk = ~np.isfinite(lons) | ~np.isfinite(lats)
    if offdisk.any():
        lons = np.where(offdisk, np.nan, lons)
        lats = np.where(offdisk, np.nan, lats)
    return lons, lats


def _ascending(source: xr.DataArray | xr.Dataset, dim: str) -> xr.DataArray | xr.Dataset:
    """xarray/scipy interpolation needs monotonically increasing coordinates."""
    values = np.asarray(source[dim].values)
    if values.size > 1 and values[0] > values[-1]:
        return source.sortby(dim)
    return source


def regrid_latlon(
    source: xr.DataArray | xr.Dataset, area: AreaDefinition
) -> xr.DataArray | xr.Dataset:
    """Bilinearly interpolate a regular lat/lon-gridded field onto ``area``.

    ``source`` must be dimensioned ``(lat, lon)`` or ``(latitude, longitude)``
    with 1-D coordinate arrays (MERRA-2 and CAMEL conventions, respectively).
    Returns the field on dims ``(y, x)`` matching the area's shape; target
    pixels outside the source extent or off the Earth disk become NaN.
    """
    lat_dim = _find_dim(source, _LAT_NAMES)
    lon_dim = _find_dim(source, _LON_NAMES)

    src = _ascending(_ascending(source, lat_dim), lon_dim)

    lons, lats = _target_lonlats(area)
    target_lat = xr.DataArray(lats, dims=("y", "x"))
    target_lon = xr.DataArray(lons, dims=("y", "x"))

    out = src.interp(
        {lat_dim: target_lat, lon_dim: target_lon},
        method="linear",
        kwargs={"bounds_error": False, "fill_value": np.nan},
    )
    return out.transpose(..., "y", "x")


def land_mask(area: AreaDefinition) -> xr.DataArray:
    """Boolean land mask (True = land) on ``area``'s grid, dims ``(y, x)``.

    Backed by the global_land_mask lookup table. Off-disk pixels are False;
    their value never matters because every BT there is NaN and the confidence
    factor propagates NaN.
    """
    lons, lats = _target_lonlats(area)
    valid = np.isfinite(lons) & np.isfinite(lats)

    mask = np.zeros(lats.shape, dtype=bool)
    if valid.any():
        # globe.is_land requires lat in [-90, 90] and lon in [-180, 180]; clip
        # away float round-off at the poles and wrap longitudes into range.
        lat_valid = np.clip(lats[valid], -90.0, 90.0)
        lon_valid = (lons[valid] + 180.0) % 360.0 - 180.0
        mask[valid] = globe.is_land(lat_valid, lon_valid)

    return xr.DataArray(mask, dims=("y", "x"))
