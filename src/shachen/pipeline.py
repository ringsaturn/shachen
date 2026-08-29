"""End-to-end DEBRA pipeline: scene + ancillary -> confidence factor fields.

Eqs. 1-22: regrid ancillary onto the scene grid, derive solar zenith and land
mask, then chain background -> cloud mask -> dust tests -> confidence.
Enhanced imagery (Eqs. 23-29) lives in :mod:`shachen.imagery`.
"""

import xarray as xr

from shachen import background as _background
from shachen import cloudmask as _cloudmask
from shachen import confidence as _confidence
from shachen import dust_tests as _dust_tests
from shachen import geo as _geo
from shachen import solar as _solar
from shachen.constants import DEFAULTS, DebraConstants

#: Variables a precomputed ``background`` Dataset must carry (the
#: :func:`shachen.background.background_signals` contract).
_BACKGROUND_VARS = (
    "rsw_bg",
    "btd_bg",
    "bt_bg_tir_86",
    "bt_bg_tir_104",
    "bt_bg_tir_123",
)

_OUTPUT_VARS = (
    "cf_comb",
    "cf_day",
    "cf_trm",
    "cf_ngt",
    "b_ngt_trm",
    "b_trm_day",
    "dt1",
    "dt2",
    "dt3",
    "cm1",
    "cm2",
    "cm3",
    "cm4",
    "cm_day",
    "cm_ngt",
    "cm_norm_day",
    "cm_norm_ngt",
    "rsw_bg",
    "btd_bg",
    "bt_bg_tir_86",
    "bt_bg_tir_104",
    "bt_bg_tir_123",
)


def _plain(field: xr.DataArray) -> xr.DataArray:
    """Strip regridding leftovers (source lat/lon coords) so merges are clean."""
    return xr.DataArray(field.data, dims=field.dims, attrs=field.attrs)


def run_debra(
    scene: xr.Dataset,
    skin_temperature: xr.DataArray,
    emissivity: xr.Dataset | None = None,
    constants: DebraConstants = DEFAULTS,
    *,
    background: xr.Dataset | None = None,
) -> xr.Dataset:
    """Run DEBRA on one scene; returns CF_comb plus all intermediate fields.

    ``scene`` is a :func:`shachen.io.satellite.load_scene` Dataset (``bt_*`` in K
    on the 2-km grid, with ``area`` and ``start_time`` attrs);
    ``skin_temperature`` is MERRA-2 TS (K) on its native lat/lon grid,
    regridded here via :func:`shachen.geo.regrid_latlon`. The visible/NIR
    reflectance variables are not used here; they feed the enhanced imagery.

    Exactly one background source must be given (ValueError otherwise):

    - ``emissivity``: the CAMEL band Dataset (``emis_*``) on its native
      lat/lon grid; regridded here, then fed through
      :func:`shachen.background.background_signals` (semianalytic mode);
    - ``background``: a precomputed Dataset **already on the scene grid**
      (e.g. :func:`shachen.composite.composite_background`) carrying
      ``rsw_bg``, ``btd_bg`` and ``bt_bg_tir_86/104/123``; missing variables
      or 2-D shapes differing from the scene raise ValueError. Its
      ``n_valid`` is passed through to the output when present.

    Returns a Dataset on the scene grid carrying ``cf_comb``, ``cf_day``,
    ``cf_trm``, ``cf_ngt``, ``cm_norm_day``, ``cm_norm_ngt``, ``dt1``-``dt3``,
    ``rsw_bg``, ``btd_bg``, and ``zenith_deg``, with the scene's ``area`` and
    ``start_time`` attrs preserved. Pixels with NaN inputs (off-disk, bad
    pixels) carry NaN confidence.
    """
    if (emissivity is None) == (background is None):
        raise ValueError(
            "exactly one background source must be given: "
            "emissivity (semianalytic) or background (precomputed)"
        )
    area = scene.attrs["area"]
    when = scene.attrs["start_time"]

    if background is not None:
        scene_shape = scene["bt_tir_104"].shape
        for name in _BACKGROUND_VARS:
            if name not in background:
                raise ValueError(f"background is missing variable {name!r}")
            if background[name].shape != scene_shape:
                raise ValueError(
                    f"shape mismatch: scene {scene_shape} "
                    f"vs background {name} {background[name].shape}"
                )

    # Ancillary onto the scene grid; drop the source lat/lon coords the
    # interpolation leaves behind so downstream merges never collide.
    ts_grid = _plain(_geo.regrid_latlon(skin_temperature, area))

    zenith = _solar.solar_zenith(area, when)
    is_land = _geo.land_mask(area)

    if background is not None:
        bg = background
    else:
        emis_grid = _geo.regrid_latlon(emissivity, area)
        emis_grid = emis_grid.map(_plain, keep_attrs=True)
        bg = _background.background_signals(ts_grid, emis_grid)
    cm = _cloudmask.cloud_mask(scene, ts_grid, constants.cloud_mask)
    dt = _dust_tests.dust_tests(scene, bg, ts_grid, is_land, constants.dust_tests)
    cf = _confidence.confidence(dt, cm, zenith, constants.confidence)

    merged = xr.merge([cf, dt, cm, bg], combine_attrs="drop")
    out = xr.Dataset({name: merged[name] for name in _OUTPUT_VARS})
    out["zenith_deg"] = zenith
    out["is_land"] = is_land
    if "n_valid" in bg:
        out["n_valid"] = bg["n_valid"]
    out.attrs.update(scene.attrs)
    return out
