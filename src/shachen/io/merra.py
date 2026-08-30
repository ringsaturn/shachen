"""MERRA-2 hourly surface skin temperature (TS from M2T1NXSLV).

Downloads via earthaccess using the Earthdata login in ``~/.netrc``. The
full tavg1_2d_slv_Nx day file is large, so it is stripped after download to a
small TS-only netCDF cache, which the pipeline regrids onto the satellite grid.
"""

import datetime as dt
from pathlib import Path

import numpy as np
import xarray as xr

SHORT_NAME = "M2T1NXSLV"

#: Single-level variables the dust-PM10 matchup keeps from M2T1NXSLV.
SURFACE_MET_SLV_VARS = ("T2M", "QV2M", "PS", "U10M", "V10M")

#: Boundary-layer height comes from the surface-flux collection.
FLX_SHORT_NAME = "M2T1NXFLX"
SURFACE_MET_FLX_VARS = ("PBLH",)


def fetch_skin_temperature(day: dt.date, out_dir: Path) -> Path:
    """Download MERRA-2 TS for ``day``; return path to a TS-only netCDF."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / f"merra2_ts_{day:%Y%m%d}.nc"
    if cache.exists():
        return cache

    import earthaccess

    earthaccess.login(strategy="netrc")
    results = earthaccess.search_data(
        short_name=SHORT_NAME,
        temporal=(day.isoformat(), day.isoformat()),
    )
    if not results:
        raise RuntimeError(f"No {SHORT_NAME} granule found for {day}")
    paths = earthaccess.download(results[:1], str(out_dir))
    full = Path(paths[0])
    with xr.open_dataset(full) as ds:
        ds[["TS"]].to_netcdf(cache)
    full.unlink()  # keep only the small TS cache
    return cache


def _download_stripped(
    short_name: str, day: dt.date, variables: tuple[str, ...], out_dir: Path
) -> xr.Dataset:
    """Download one MERRA-2 day granule and return only ``variables``."""
    import earthaccess

    earthaccess.login(strategy="netrc")
    results = earthaccess.search_data(
        short_name=short_name,
        temporal=(day.isoformat(), day.isoformat()),
    )
    if not results:
        raise RuntimeError(f"No {short_name} granule found for {day}")
    paths = earthaccess.download(results[:1], str(out_dir))
    full = Path(paths[0])
    with xr.open_dataset(full) as ds:
        stripped = ds[list(variables)].load()
    full.unlink()  # keep only the small stripped cache
    return stripped


def fetch_surface_meteorology(day: dt.date, out_dir: Path) -> Path:
    """Download the dust-PM10 met covariates for ``day``; return a cache path.

    The cache ``merra2_met_<%Y%m%d>.nc`` merges ``SURFACE_MET_SLV_VARS`` from
    M2T1NXSLV with ``SURFACE_MET_FLX_VARS`` (PBLH) from M2T1NXFLX, both
    hourly on the native 0.5 x 0.625 degree grid.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / f"merra2_met_{day:%Y%m%d}.nc"
    if cache.exists():
        return cache

    slv = _download_stripped(SHORT_NAME, day, SURFACE_MET_SLV_VARS, out_dir)
    flx = _download_stripped(FLX_SHORT_NAME, day, SURFACE_MET_FLX_VARS, out_dir)
    xr.merge([slv, flx]).to_netcdf(cache)
    return cache


def _interp_to_time(dataset: xr.Dataset, when: dt.datetime) -> xr.Dataset:
    """Linear interpolation in time, clamped to the file's own time axis.

    MERRA-2 hourly means are stamped at half past the hour (00:30, 01:30,
    ...), so a scan on the whole hour at either end of a day file — 00:00
    UTC, say — falls outside that axis. Plain ``interp`` extrapolates to
    NaN there, and a NaN skin temperature silently empties every field
    downstream of the background and the cloud mask. Clamping reads the
    nearest hourly mean instead (00:00 gets the 00:30 mean, half an hour
    off), which is the same approximation the interpolation already makes
    between stamps, and never invents a value the day never had.
    """
    stamp = np.datetime64(when.replace(tzinfo=None), "ns")
    times = dataset["time"].values
    return dataset.interp(time=min(max(stamp, times.min()), times.max()))


def load_surface_meteorology(path: Path, when: dt.datetime) -> xr.Dataset:
    """Load the met covariates interpolated to ``when``.

    Same half-past-the-hour time stamps, and the same clamping at the ends
    of the day, as :func:`load_skin_temperature`.
    """
    with xr.open_dataset(path) as ds:
        return _interp_to_time(ds, when).load()


def load_skin_temperature(path: Path, when: dt.datetime) -> xr.DataArray:
    """Load TS (K) interpolated to ``when`` on the MERRA-2 grid.

    Linear in time between the half-past-the-hour stamps, clamped to the
    file's first and last hourly mean — see :func:`_interp_to_time`.
    """
    with xr.open_dataset(path) as ds:
        ts = _interp_to_time(ds[["TS"]], when).load()["TS"]
    ts.attrs["units"] = "K"
    return ts
