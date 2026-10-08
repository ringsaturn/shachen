"""Satellite L1b -> calibrated BT/reflectance fields on the 2-km IR grid.

satpy does the radiance -> BT / reflectance calibration and the native
down-sampling of the VIS/NIR bands onto the coarsest (2 km) grid, so ABI
netCDF and Himawari AHI HSD share one code path.
"""

import tempfile
from collections.abc import Iterable
from pathlib import Path

import satpy
import xarray as xr
from satpy import Scene

from shachen.constants import ABI_BANDS, AHI_BANDS, Band

#: Bands calibrated to reflectance (%); the rest to brightness temperature (K).
_REFLECTIVE = {Band.VIS_064, Band.NIR_160}

_READERS = {"abi_l1b": ABI_BANDS, "ahi_hsd": AHI_BANDS}


def load_scene(
    files: Iterable[str | Path],
    reader: str = "abi_l1b",
    roles: Iterable[Band] | None = None,
    bbox: tuple[float, float, float, float] | None = None,
) -> xr.Dataset:
    """Load L1b files and return a Dataset on the sensor's 2-km fixed grid.

    Variables are named by DEBRA role: ``bt_tir_104`` etc. for emissive
    bands (K), ``refl_vis_064`` etc. for reflective bands (%). The
    pyresample AreaDefinition is attached as ``ds.attrs["area"]`` and the
    scan start time as ``ds.attrs["start_time"]``.

    ``roles`` restricts loading to a subset of DEBRA bands (e.g.
    ``shachen.composite.COMPOSITE_BANDS`` for TIR-only composite days, whose
    directories hold only those bands' files); ``None`` loads all seven.

    ``bbox`` = ``(lon_min, lat_min, lon_max, lat_max)`` in degrees crops the
    scene (satpy ``Scene.crop(ll_bbox=...)``) *after* the native resample,
    so every variable and ``attrs["area"]`` describe the cropped grid —
    This trims AHI full disks (5500 x 5500 at 2 km) to the domain of interest.
    ``None`` keeps the full scene.
    """
    band_map = _READERS[reader]
    if roles is not None:
        band_map = {role: band_map[role] for role in roles}
    # satpy unpacks bz2-compressed AHI HSD segments into ``config["tmp_dir"]``
    # and relies on ``weakref.finalize(self, self._cleanup)`` to delete them --
    # a bound method that keeps the file handler alive, so the files outlive
    # the Scene and are only removed at interpreter exit (never, if the
    # process is killed). Measured 2026-09-17: 3 leaked files per single-band
    # load, ~430 MB per 15-scene composite, 93 GB on one batch worker. Giving
    # each load its own tmp_dir and computing the arrays before it goes away
    # bounds that to the life of this call; ABI netCDF never touches it.
    with (
        tempfile.TemporaryDirectory(prefix="shachen-l1b-") as tmp_dir,
        satpy.config.set(tmp_dir=tmp_dir),
    ):
        scn = Scene(filenames=[str(f) for f in files], reader=reader)
        scn.load(list(band_map.values()))
        # Native resampling aggregates the 0.5/1 km reflective bands onto the
        # coarsest-loaded (2 km IR) grid without interpolation artifacts.
        scn = scn.resample(scn.coarsest_area(), resampler="native")
        if bbox is not None:
            # Cropping after the resample keeps every dataset on one shared area.
            scn = scn.crop(ll_bbox=bbox)

        data_vars = {}
        for role, name in band_map.items():
            da = scn[name].drop_vars("crs", errors="ignore").load()
            prefix = "refl" if role in _REFLECTIVE else "bt"
            data_vars[f"{prefix}_{role.value}"] = da
    ds = xr.Dataset(data_vars)
    ds.attrs["area"] = scn.coarsest_area()
    ds.attrs["start_time"] = scn.start_time
    ds.attrs["reader"] = reader
    return ds
