"""End-to-end entry points: one call per algorithm, scene in, fields out.

:func:`run_debra` is DEBRA, Eqs. 1-22: regrid ancillary onto the scene grid,
derive solar zenith and land mask, then chain background -> cloud mask ->
dust tests -> confidence. Enhanced imagery (Eqs. 23-29) is in
:mod:`shachen.imagery`.

:func:`run_zhouye` is ZHOUYE, the diurnally consistent confidence scheme
built on DEBRA's tests and cloud mask: the same chain with the
:data:`shachen.constants.ZHOUYE` preset, tagged ``scheme="zhouye"``.

:func:`run_dust_rgb` is the classic Dust RGB baseline, which needs no
ancillary data at all but does need to know which sensor it is looking at.
Callers reach all of them through this module rather than the per-equation
modules, so adding an algorithm does not change how one is run.
"""

import xarray as xr

from shachen import background as _background
from shachen import cloudmask as _cloudmask
from shachen import confidence as _confidence
from shachen import dust_tests as _dust_tests
from shachen import dustrgb as _dustrgb
from shachen import geo as _geo
from shachen import solar as _solar
from shachen.constants import (
    DEFAULTS,
    DUST_RGB,
    DUST_RGB_BY_READER,
    ZHOUYE,
    DebraConstants,
    DustRGBConstants,
)

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
    if "dt2_fixed" in dt:
        # ZHOUYE: the second DT2 reading that Eqs. 17-18 used; dt2 stays Eq. 14.
        out["dt2_fixed"] = dt["dt2_fixed"]
    out["zenith_deg"] = zenith
    out["is_land"] = is_land
    if "n_valid" in bg:
        out["n_valid"] = bg["n_valid"]
    out.attrs.update(scene.attrs)
    # Which confidence scheme produced cf_comb. DEBRA unless the constants
    # turn on a ZHOUYE hook; run_zhouye sets it explicitly.
    out.attrs["scheme"] = (
        "zhouye"
        if constants.dust_tests.dt2_fixed_interval is not None
        or constants.confidence.corroborate_dt3
        else "debra"
    )
    return out


def run_zhouye(
    scene: xr.Dataset,
    skin_temperature: xr.DataArray,
    emissivity: xr.Dataset | None = None,
    constants: DebraConstants = ZHOUYE,
    *,
    background: xr.Dataset | None = None,
) -> xr.Dataset:
    """Run ZHOUYE on one scene: :func:`run_debra`'s chain with the ZHOUYE preset.

    Same inputs and outputs as :func:`run_debra` (it *is* that chain -- the
    tests, cloud mask, blending and imagery are DEBRA's), with three things
    switched on through the constants: DT2 on a fixed interval in Eqs. 17-18
    instead of Eq. 14, DT3 needing corroboration there, and terminator and
    night Eq. 19 intervals fitted so the same dust reads the same as by day.
    Eq. 16 is untouched, so by day (solar zenith below 75 deg) ``cf_comb``
    equals a :data:`ABI_TUNED` run's bit for bit; the two diverge across the
    terminator and at night. See
    :data:`shachen.constants.ZHOUYE` for the claim each hook buys and what
    it does not. The output is tagged ``attrs["scheme"] == "zhouye"`` and
    carries ``dt2_fixed`` next to ``dt2``.
    """
    out = run_debra(scene, skin_temperature, emissivity, constants, background=background)
    out.attrs["scheme"] = "zhouye"
    return out


def run_dust_rgb(scene: xr.Dataset, constants: DustRGBConstants | None = None) -> xr.Dataset:
    """Run the classic Dust RGB baseline on one scene.

    The counterpart of :func:`run_debra` for the recipe in
    :mod:`shachen.dustrgb`: same scene in, but no ancillary data, no cloud mask
    and no confidence field — three fixed stretches of ``bt_tir_86/104/112/123``
    (11.2 um is the extra band DEBRA itself never reads). A scene loaded with
    ``roles=DEBRA_BANDS`` therefore raises ValueError here.

    The stretches are per sensor. Unlike DEBRA, the Dust RGB has no one
    canonical set of numbers: it was tuned for SEVIRI and then re-tuned for each
    later imager, because the corresponding channels do not sit at the same
    wavelengths. With ``constants=None`` (the default) the set is chosen from
    ``scene.attrs["reader"]`` through :data:`shachen.constants.DUST_RGB_BY_READER`
    — ABI gets the Quick Guide's adjusted values, AHI the original SEVIRI ones —
    so the baseline matches that sensor's operational product. An unknown or
    absent reader falls
    back to :data:`shachen.constants.DUST_RGB` (SEVIRI); pass ``constants``
    explicitly to pin one set across sensors, e.g. to compare the two.

    Returns a Dataset carrying ``dust_rgb`` — dims ``(y, x, gun)``, floats in
    [0, 1], ready for :func:`shachen.imagery.to_uint8` — with the scene's
    ``area`` and ``start_time`` attrs preserved, so it merges straight into a
    :func:`run_debra` result for side-by-side rendering.
    """
    if constants is None:
        constants = DUST_RGB_BY_READER.get(scene.attrs.get("reader"), DUST_RGB)
    out = xr.Dataset({"dust_rgb": _dustrgb.dust_rgb(scene, constants)})
    out.attrs.update(scene.attrs)
    return out
