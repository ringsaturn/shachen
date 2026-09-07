"""DEBRA confidence factor, Eqs. 16-22 (Eqs. 21-22 per the 2020 erratum).

Combines the dust tests under cloud-mask suppression into day / terminator /
night confidence factors (Eqs. 16-18), normalizes each (Eq. 19), and blends
them across the terminator with solar-zenith weights ``B_ngt_trm`` and
``B_trm_day`` (Eqs. 20-21, exponent 1.5) into the final ``CF_comb`` in [0, 1]
(Eq. 22).

Eq. 19 as printed uses one interval, (0.25, 2.50), for all three branches.
That is a scale error, and a visible one: Eq. 16 sums three DT terms and can
reach 3.0, Eq. 18 sums two and can reach 1.5, so one shared interval caps
CF_ngt at (1.5 - 0.25) / (2.50 - 0.25) = 0.556 while CF_day reaches 1. Dust
appears to fade at dusk. The deviation here is to give each branch its own
interval (:class:`shachen.constants.ConfidenceConstants`), with the
terminator interpolating between them on the Eq. 20 weight.
"""

import numpy as np
import xarray as xr

from shachen.constants import DEFAULTS, ConfidenceConstants
from shachen.norm import normalize, normalize_cos_zenith, normalize_interp

_TEST_VARS = ("dt1", "dt2", "dt3")
_CLOUD_VARS = ("cm_norm_day", "cm_norm_ngt")


def _require(dataset: xr.Dataset, names: tuple[str, ...], label: str) -> None:
    missing = [name for name in names if name not in dataset]
    if missing:
        raise ValueError(f"{label} is missing required variable(s): {missing}")


def _check_shapes(fields: dict[str, xr.DataArray]) -> None:
    shapes = {name: tuple(field.shape) for name, field in fields.items()}
    if len(set(shapes.values())) > 1:
        raise ValueError(f"all 2-D inputs must share one shape, got {shapes}")


def confidence(
    tests: xr.Dataset,
    cloud: xr.Dataset,
    zenith_deg: xr.DataArray,
    constants: ConfidenceConstants = DEFAULTS.confidence,
) -> xr.Dataset:
    """Confidence factors and blend weights on the scene grid.

    ``tests`` needs ``dt1``, ``dt2``, ``dt3`` (from
    :func:`shachen.dust_tests.dust_tests`); ``cloud`` needs ``cm_norm_day``,
    ``cm_norm_ngt`` (from :func:`shachen.cloudmask.cloud_mask`); ``zenith_deg``
    is the solar zenith angle in degrees. All 2-D inputs must share one shape
    (ValueError otherwise); NaN propagates.

    Returns a Dataset with ``cf_day``, ``cf_trm``, ``cf_ngt`` (each already
    normalized per Eq. 19; CF_trm and CF_ngt use ``cm_norm_day`` and
    ``cm_norm_ngt`` respectively, CF_ngt takes ``max(DT1, DT2)``), the blend
    weights ``b_ngt_trm``, ``b_trm_day`` (via
    :func:`shachen.norm.normalize_cos_zenith`), and ``cf_comb`` (Eq. 22).

    The Eq. 19 intervals are per branch: ``constants.cf_norm_day`` for
    CF_day, ``constants.cf_norm_ngt`` for CF_ngt, and for CF_trm the two
    interpolated on ``b_ngt_trm``, so its scale matches whichever neighbour
    Eq. 22 is blending it into. Passing ``cf_norm=`` to
    :class:`shachen.constants.ConfidenceConstants` sets both intervals to one
    value and restores the printed single-interval behaviour.
    """
    _require(tests, _TEST_VARS, "tests")
    _require(cloud, _CLOUD_VARS, "cloud")

    dt1 = tests["dt1"]
    dt2 = tests["dt2"]
    dt3 = tests["dt3"]
    cm_day = cloud["cm_norm_day"]
    cm_ngt = cloud["cm_norm_ngt"]

    _check_shapes(
        {
            "dt1": dt1,
            "dt2": dt2,
            "dt3": dt3,
            "cm_norm_day": cm_day,
            "cm_norm_ngt": cm_ngt,
            "zenith_deg": zenith_deg,
        }
    )

    c = constants

    # Eq. 16: full daytime confidence, all three dust tests at full weight.
    cf_day_raw = (dt1 + dt2 + dt3) * (1.0 - cm_day)
    # Eq. 17: terminator confidence, DT3 down-weighted, daytime cloud mask.
    cf_trm_raw = (dt1 + dt2 + c.dt3_weight_trm * dt3) * (1.0 - cm_day)
    # Eq. 18: night confidence, max(DT1, DT2), nighttime cloud mask.
    cf_ngt_raw = (np.maximum(dt1, dt2) + c.dt3_weight_ngt * dt3) * (1.0 - cm_ngt)

    # Eqs. 20-21 (21 per the 2020 erratum): cos-zenith blend weights.
    b_ngt_trm = normalize_cos_zenith(zenith_deg, c.ngt_trm_zenith_deg, c.blend_exponent)
    b_trm_day = normalize_cos_zenith(zenith_deg, c.trm_day_zenith_deg, c.blend_exponent)

    # Eq. 19: normalize each variant onto [0, 1]. Day and night use their own
    # intervals because Eqs. 16 and 18 have different raw ceilings (3.0 and
    # 1.5); the terminator, whose ceiling is between them, rides the Eq. 20
    # weight from one interval to the other.
    cf_day = normalize(cf_day_raw, c.cf_norm_day)
    cf_trm = normalize_interp(cf_trm_raw, c.cf_norm_day, c.cf_norm_ngt, b_ngt_trm)
    cf_ngt = normalize(cf_ngt_raw, c.cf_norm_ngt)

    # Eq. 22 (erratum): nested day / terminator / night blend.
    cf_comb = b_trm_day * cf_day + (1.0 - b_trm_day) * (
        b_ngt_trm * cf_trm + (1.0 - b_ngt_trm) * cf_ngt
    )

    return xr.Dataset(
        {
            "cf_day": cf_day,
            "cf_trm": cf_trm,
            "cf_ngt": cf_ngt,
            "b_ngt_trm": b_ngt_trm,
            "b_trm_day": b_trm_day,
            "cf_comb": cf_comb,
        }
    )
