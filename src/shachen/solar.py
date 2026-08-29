"""Per-pixel solar zenith angle for the day/terminator/night blends.

Feeds the confidence-factor blending weights (Eqs. 20-21) and the
baseline-image blend (Eq. 25). Implemented with pyorbital's astronomy
routines on the area's lon/lat arrays.
"""

import datetime as dt

import numpy as np
import xarray as xr
from pyorbital import astronomy
from pyresample.geometry import AreaDefinition


def solar_zenith(area: AreaDefinition, when: dt.datetime) -> xr.DataArray:
    """Solar zenith angle in degrees on ``area``'s grid, dims ``(y, x)``.

    ``when`` is the (naive UTC) observation time; the scene start time is
    accurate enough for the 30-deg-wide terminator blend. Off-disk pixels
    (non-finite lon/lat) are NaN.
    """
    lons, lat = area.get_lonlats()
    lons = np.asarray(lons, dtype=np.float64)
    lats = np.asarray(lat, dtype=np.float64)

    # Off-disk pixels come back as +/-inf (or already NaN); make them NaN so the
    # arccos below propagates NaN instead of raising or clipping to a bogus angle.
    offdisk = ~(np.isfinite(lons) & np.isfinite(lats))
    lons = np.where(offdisk, np.nan, lons)
    lats = np.where(offdisk, np.nan, lats)

    with np.errstate(invalid="ignore"):
        cos_zen = astronomy.cos_zen(when, lons, lats)
        zenith = np.degrees(np.arccos(np.clip(cos_zen, -1.0, 1.0)))
    zenith = np.where(offdisk, np.nan, zenith)

    return xr.DataArray(zenith, dims=("y", "x"), name="solar_zenith_angle")
