"""DEBRA enhanced imagery, Eqs. 23-29 (Eqs. 24-25 per the 2020 erratum).

Builds the day/night blended baseline image (Eqs. 23-26) and modulates it by
the combined confidence factor in each color gun (Eqs. 27-29). Pure array
numerics; map rendering is in :mod:`shachen.render`.
"""

import numpy as np
import xarray as xr

from shachen.constants import (
    COLOR_DIMMING,
    DEFAULTS,
    Bounds,
    DebraConstants,
    ImageryConstants,
)
from shachen.norm import normalize

_SCENE_VARS = ("refl_vis_064", "bt_tir_104")
_DEBRA_VARS = ("cf_comb", "zenith_deg")

#: Color guns of the Eqs. 27-29 composite, in the coordinate's order.
_GUNS = ("r", "g", "b")


def _require(dataset: xr.Dataset, names: tuple[str, ...], label: str) -> None:
    missing = [name for name in names if name not in dataset]
    if missing:
        raise ValueError(f"{label} is missing required variable(s): {missing}")


def _check_shapes(fields: dict[str, xr.DataArray]) -> None:
    shapes = {name: tuple(field.shape) for name, field in fields.items()}
    if len(set(shapes.values())) > 1:
        raise ValueError(f"all 2-D inputs must share one shape, got {shapes}")


def _domain_normalized(field: xr.DataArray, invert: bool = False) -> xr.DataArray:
    """Domain min/max normalization of ``field`` (Eq. 23; Eq. 24 if inverted).

    The min/max are taken over the whole domain skipping NaN. A degenerate
    domain (no finite pixel, or max == min) yields 0.0 everywhere rather than
    a division by zero.
    """
    values = np.asarray(field.values, dtype=float)
    finite = np.isfinite(values)
    if finite.any():
        lo = values[finite].min()
        hi = values[finite].max()
    else:
        lo = hi = 0.0

    if hi == lo:
        # Degenerate domain: the component carries no information.
        scaled = np.zeros_like(values)
    else:
        scaled = (values - lo) / (hi - lo)
        if invert:
            scaled = 1.0 - scaled

    return xr.DataArray(scaled, dims=field.dims, coords=field.coords)


def baseline_image(
    scene: xr.Dataset,
    zenith_deg: xr.DataArray,
    constants: ImageryConstants = DEFAULTS.imagery,
) -> xr.Dataset:
    """Day/night blended baseline image BI (Eqs. 23-26).

    ``scene`` needs ``refl_vis_064`` (any linear reflectance unit; the Eq. 23
    domain min/max normalization cancels units) and ``bt_tir_104`` (K);
    ``zenith_deg`` is the solar zenith angle in degrees. Missing variables and
    2-D shape mismatches raise ValueError.

    Returns a Dataset on dims ``(y, x)`` with

    - ``vis_bg`` (Eq. 23): VIS normalized by the domain min/max reflectance;
    - ``ir_bg`` (Eq. 24, erratum): inverted normalized BT10.4, 1.0 at the
      domain cold bound;
    - ``b_bg`` (Eq. 25, erratum): ``1 - N(zenith; 79, 89)**1.5`` evaluated in
      zenith-*degree* space (unlike the cos-space Eqs. 20-21);
    - ``bi`` (Eq. 26): ``b_bg*vis_bg + (1 - b_bg)*ir_bg``, in [0, 1].

    NaN semantics: NaN pixels propagate, except a component whose Eq. 26
    blend weight is exactly zero is ignored (a NaN VIS pixel on the night
    side leaves ``bi = ir_bg``). A degenerate component domain (all-NaN, or
    domain max == min, e.g. an all-dark night VIS band) yields 0.0 for that
    component everywhere instead of dividing by zero.
    """
    _require(scene, _SCENE_VARS, "scene")

    refl = scene["refl_vis_064"]
    bt = scene["bt_tir_104"]
    _check_shapes({"refl_vis_064": refl, "bt_tir_104": bt, "zenith_deg": zenith_deg})

    c = constants

    # Eq. 23: VIS normalized by the domain reflectance range.
    vis_bg = _domain_normalized(refl)
    # Eq. 24 (erratum): inverted normalized BT10.4, 1.0 at the domain cold bound.
    ir_bg = _domain_normalized(bt, invert=True)
    # Eq. 25 (erratum): blend weight in zenith-*degree* space (not cos space).
    b_bg = 1.0 - normalize(zenith_deg, c.bg_blend_zenith_deg) ** c.bg_blend_exponent

    # Eq. 26, with zero-weight components excluded so that a NaN pixel on the
    # unused side of the terminator does not poison the blend.
    weight = np.asarray(b_bg.values, dtype=float)
    vis = np.asarray(vis_bg.values, dtype=float)
    ir = np.asarray(ir_bg.values, dtype=float)
    blended = np.where(
        weight == 0.0,
        ir,
        np.where(weight == 1.0, vis, weight * vis + (1.0 - weight) * ir),
    )
    bi = xr.DataArray(np.clip(blended, 0.0, 1.0), dims=vis_bg.dims, coords=vis_bg.coords)

    return xr.Dataset({"vis_bg": vis_bg, "ir_bg": ir_bg, "b_bg": b_bg, "bi": bi})


def enhanced_rgb(
    baseline: xr.DataArray,
    cf_comb: xr.DataArray,
    constants: ImageryConstants = DEFAULTS.imagery,
    dimming: tuple[float, float, float] | None = None,
) -> xr.DataArray:
    """Color-modulated composite (Eqs. 27-29 + the Eq. 3 gun rescale).

    Each gun is ``BI*(1 - min(CF, cf_cap)) + D_gun*CF`` with the per-gun
    dimming triple ``dimming`` (``None`` selects ``COLOR_DIMMING["yellow"]``,
    i.e. D = 1.0 on red/green and ``blue_dimming`` on blue), then rescaled
    with ``N(gun; 0, gun_max)`` onto [0, 1].

    Returns a float DataArray with dims ``("y", "x", "gun")`` and coordinate
    ``gun = ["r", "g", "b"]`` (the dim is not named ``rgb`` so the array can
    live in a Dataset as variable ``rgb`` without a name collision). NaN
    propagates; a shape mismatch between ``baseline`` and ``cf_comb`` raises
    ValueError.
    """
    _check_shapes({"baseline": baseline, "cf_comb": cf_comb})

    c = constants
    if dimming is None:
        dimming = COLOR_DIMMING["yellow"]

    bi = np.asarray(baseline.values, dtype=float)
    cf = np.asarray(cf_comb.values, dtype=float)

    # Eqs. 27-29: the confidence factor bleeds into each gun, capped at cf_cap
    # in the baseline term so dust never fully erases the scene.
    damped = bi * (1.0 - np.minimum(cf, c.cf_cap))
    gun_bounds = Bounds(0.0, c.gun_max)
    guns = [normalize(damped + gun_dimming * cf, gun_bounds) for gun_dimming in dimming]

    return xr.DataArray(
        np.stack(guns, axis=-1),
        dims=tuple(baseline.dims) + ("gun",),
        coords={**dict(baseline.coords), "gun": list(_GUNS)},
    )


def to_uint8(rgb: xr.DataArray) -> np.ndarray:
    """Convert the [0, 1] float RGB composite to a ``(y, x, 3)`` uint8 array.

    ``rgb`` is an :func:`enhanced_rgb`-shaped DataArray (gun axis last).
    Values are scaled by 255 and rounded with :func:`numpy.rint`; NaN pixels
    (off-disk) become 0 (black).
    """
    values = np.rint(np.asarray(rgb.values, dtype=float) * 255.0)
    values = np.where(np.isfinite(values), values, 0.0)
    return np.clip(values, 0.0, 255.0).astype(np.uint8)


def debra_imagery(
    scene: xr.Dataset,
    debra: xr.Dataset,
    constants: DebraConstants = DEFAULTS,
    dimming: tuple[float, float, float] | None = None,
) -> xr.Dataset:
    """Full imagery chain (Eqs. 23-29) from a scene and a run_debra output.

    ``scene`` needs ``refl_vis_064`` and ``bt_tir_104``; ``debra`` needs
    ``cf_comb`` and ``zenith_deg`` (from :func:`shachen.pipeline.run_debra`).
    Missing variables raise ValueError.

    Returns a Dataset with the :func:`baseline_image` fields (``vis_bg``,
    ``ir_bg``, ``b_bg``, ``bi``) plus ``rgb`` from :func:`enhanced_rgb`,
    with ``scene``'s attrs preserved.
    """
    _require(debra, _DEBRA_VARS, "debra")

    out = baseline_image(scene, debra["zenith_deg"], constants.imagery)
    out["rgb"] = enhanced_rgb(out["bi"], debra["cf_comb"], constants.imagery, dimming)
    out.attrs.update(scene.attrs)
    return out
