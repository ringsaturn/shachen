"""IR surface emissivity climatology, interpolated to sensor band centers.

The paper uses the monthly UW Baseline Fit (UWBF, CIMSS) 0.05-degree
climatology, which needs a separate CIMSS registration. The default here is
its successor, CAMEL (CAM5K30EM on NASA LP DAAC), reachable with the same
Earthdata credentials as MERRA-2. A downloaded UWBF file works too: both
formats carry an emissivity cube on labelled hinge-point wavelengths, and
interpolation is linear in wavelength, as the paper implies.
"""

import datetime as dt
import re
from pathlib import Path

import numpy as np
import xarray as xr

from shachen.constants import BAND_CENTER_UM, Band

SHORT_NAME = "CAM5K30EM"  # CAMEL monthly 0.05-deg emissivity

#: DEBRA only needs emissivity for the emissive bands used by the background.
EMISSIVE_BANDS = (Band.SWIR_39, Band.WV_62, Band.TIR_86, Band.TIR_104, Band.TIR_123)


def fetch_emissivity(month: dt.date, out_dir: Path) -> Path:
    """Download the CAMEL monthly emissivity file covering ``month``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    first = month.replace(day=1)
    existing = sorted(out_dir.glob(f"CAM5K30EM*{first:%Y%m}*.nc"))
    if existing:
        return existing[0]

    import earthaccess

    earthaccess.login(strategy="netrc")
    # Search the whole month, then pick the granule whose native filename
    # carries this month's YYYYMM token: a point-in-time temporal search can
    # match the *previous* month's granule, whose coverage interval ends on
    # the 1st (seen with 2017-03: the Feb file spans Feb 1 - Mar 1).
    last = (first.replace(day=28) + dt.timedelta(days=4)).replace(day=1) - dt.timedelta(days=1)
    results = earthaccess.search_data(
        short_name=SHORT_NAME,
        temporal=(first.isoformat(), last.isoformat()),
    )
    token = f"{first:%Y%m}"
    matched = [g for g in results if token in " ".join(g.data_links())]
    if not matched:
        raise RuntimeError(
            f"No {SHORT_NAME} granule matching {token} found "
            f"({len(results)} granules in the temporal window)"
        )
    paths = earthaccess.download(matched[:1], str(out_dir))
    path = Path(paths[0])
    # Normalize the name so the cache glob above finds it next time.
    canonical = out_dir / f"{SHORT_NAME}_{first:%Y%m}.nc"
    if path != canonical:
        path.rename(canonical)
    return canonical


def _find_emissivity_cube(ds: xr.Dataset) -> tuple[xr.DataArray, np.ndarray, str]:
    """Return (emissivity[cube], hinge_wavelengths_um, spectral_dim_name)."""
    for var in ("camel_emis", "camel_emissivity", "CAMEL_emissivity", "emis", "emissivity"):
        if var in ds:
            da = ds[var]
            break
    else:
        raise KeyError(f"No emissivity variable found; file has {list(ds.data_vars)}")
    spectral_dims = [d for d in da.dims if d not in ("latitude", "longitude", "lat", "lon")]
    if len(spectral_dims) != 1:
        raise ValueError(f"Cannot identify spectral dim among {da.dims}")
    dim = spectral_dims[0]
    for wl_var in ("wavelength", "wavelengths", dim):
        if wl_var in ds.variables:
            wl = np.asarray(ds[wl_var].values, dtype=float)
            break
    else:
        # CAM5K30EM carries the hinge points only in the variable's comment
        # string ("Emissivity at 3.6, 4.3, ... 14.3 micron"), not as a coord.
        comment = str(da.attrs.get("comment", ""))
        wl = np.array([float(m) for m in re.findall(r"\d+\.\d+", comment)])
        if wl.size != da.sizes[dim]:
            raise KeyError(
                f"No hinge-point wavelength coordinate found and the comment "
                f"attr yields {wl.size} values for spectral dim of size "
                f"{da.sizes[dim]}"
            )
    return da, wl, dim


def load_band_emissivity(path: Path, bands: tuple[Band, ...] = EMISSIVE_BANDS) -> xr.Dataset:
    """Load hinge-point emissivity and interpolate to DEBRA band centers.

    Returns a Dataset with one ``emis_<band>`` variable per requested band
    on the climatology's native lat/lon grid, values masked where the file
    has no retrieval (ocean/fill).
    """
    with xr.open_dataset(path) as ds:
        da, wl, dim = _find_emissivity_cube(ds)
        da = da.load()

    scale = da.attrs.get("scale_factor", 1.0)
    fill = da.attrs.get("_FillValue", None)
    values = da.values.astype("float32")
    if fill is not None:
        values = np.where(values == fill, np.nan, values)
    values = values * scale
    da = da.copy(data=values)

    out = {}
    for band in bands:
        center = BAND_CENTER_UM[band]
        # Linear in wavelength between hinge points; clamp outside the hinges.
        center = float(np.clip(center, wl.min(), wl.max()))
        idx = np.searchsorted(wl, center)
        if idx == 0 or wl[idx - 1] == center:
            e = da.isel({dim: max(idx - 1, 0)})
        else:
            lo, hi = wl[idx - 1], wl[idx]
            w = (center - lo) / (hi - lo)
            e = da.isel({dim: idx - 1}) * (1 - w) + da.isel({dim: idx}) * w
        e = e.clip(0.0, 1.0)
        e.attrs = {"band_center_um": BAND_CENTER_UM[band]}
        out[f"emis_{band.value}"] = e
    return xr.Dataset(out)
