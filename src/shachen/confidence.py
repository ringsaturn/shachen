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

The two halves are also callable on their own: :func:`confidence_raw` is
Eqs. 16-18, :func:`confidence_norm` is Eqs. 19-22. Everything the intervals
touch lives in the second half, so a stored raw sum can be re-normalized on a
different interval without the dust tests or the L1b behind them.
"""

import numpy as np
import xarray as xr

from shachen.constants import DEFAULTS, ConfidenceConstants
from shachen.norm import normalize, normalize_cos_zenith, normalize_interp

_TEST_VARS = ("dt1", "dt2", "dt3")
_CLOUD_VARS = ("cm_norm_day", "cm_norm_ngt")
_RAW_VARS = ("cf_day_raw", "cf_trm_raw", "cf_ngt_raw")


def _require(dataset: xr.Dataset, names: tuple[str, ...], label: str) -> None:
    missing = [name for name in names if name not in dataset]
    if missing:
        raise ValueError(f"{label} is missing required variable(s): {missing}")


def _check_shapes(fields: dict[str, xr.DataArray]) -> None:
    shapes = {name: tuple(field.shape) for name, field in fields.items()}
    if len(set(shapes.values())) > 1:
        raise ValueError(f"all 2-D inputs must share one shape, got {shapes}")


def confidence_raw(
    tests: xr.Dataset,
    cloud: xr.Dataset,
    constants: ConfidenceConstants = DEFAULTS.confidence,
) -> xr.Dataset:
    """Eqs. 16-18 before the Eq. 19 normalization.

    Returns ``cf_day_raw``, ``cf_trm_raw`` and ``cf_ngt_raw`` on the scene
    grid, in the units the DT tests come in (each already suppressed by its
    cloud mask). Eq. 19 is a clipped monotone map of these, so they are what
    an evaluation needs in order to re-normalize with different intervals
    later without re-running the retrieval -- and, because a monotone map
    commutes with taking a median, they can be spatially aggregated first
    and normalized afterwards, which the normalized factors cannot.

    ``tests`` needs ``dt1``, ``dt2``, ``dt3``; ``cloud`` needs
    ``cm_norm_day``, ``cm_norm_ngt``.
    """
    _require(tests, _TEST_VARS, "tests")
    _require(cloud, _CLOUD_VARS, "cloud")

    dt1 = tests["dt1"]
    dt2 = tests["dt2"]
    dt3 = tests["dt3"]
    cm_day = cloud["cm_norm_day"]
    cm_ngt = cloud["cm_norm_ngt"]
    c = constants

    return xr.Dataset(
        {
            # Eq. 16: full daytime confidence, all three tests at full weight.
            "cf_day_raw": (dt1 + dt2 + dt3) * (1.0 - cm_day),
            # Eq. 17: terminator, DT3 down-weighted, daytime cloud mask.
            "cf_trm_raw": (dt1 + dt2 + c.dt3_weight_trm * dt3) * (1.0 - cm_day),
            # Eq. 18: night, max(DT1, DT2), nighttime cloud mask.
            "cf_ngt_raw": (np.maximum(dt1, dt2) + c.dt3_weight_ngt * dt3) * (1.0 - cm_ngt),
        }
    )


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

    return confidence_norm(confidence_raw(tests, cloud, constants), zenith_deg, constants)


def blend_confidence(cf_day, cf_trm, cf_ngt, b_ngt_trm, b_trm_day):
    """Eq. 22 (erratum): the nested day / terminator / night blend.

    Split out so that anything normalizing the branches differently --
    :mod:`shachen.calibration` does -- blends them the same way rather than
    keeping a second copy of Eq. 22.
    """
    return b_trm_day * cf_day + (1.0 - b_trm_day) * (
        b_ngt_trm * cf_trm + (1.0 - b_ngt_trm) * cf_ngt
    )


def confidence_norm(
    raw: xr.Dataset,
    zenith_deg,
    constants: ConfidenceConstants = DEFAULTS.confidence,
) -> xr.Dataset:
    """Eqs. 19-22 on the raw sums :func:`confidence_raw` produces.

    ``raw`` needs ``cf_day_raw``, ``cf_trm_raw``, ``cf_ngt_raw``. Split out
    from :func:`confidence` because everything the Eq. 19 intervals touch is
    here and nothing upstream of them is: re-tuning the intervals on stored
    raw sums is this call, and needs neither the L1b nor the dust tests.

    Returns the same variables :func:`confidence` does.
    """
    _require(raw, _RAW_VARS, "raw")
    c = constants

    # Eqs. 20-21 (21 per the 2020 erratum): cos-zenith blend weights.
    b_ngt_trm = normalize_cos_zenith(zenith_deg, c.ngt_trm_zenith_deg, c.blend_exponent)
    b_trm_day = normalize_cos_zenith(zenith_deg, c.trm_day_zenith_deg, c.blend_exponent)

    # Eq. 19: normalize each variant onto [0, 1]. Day and night use their own
    # intervals because Eqs. 16 and 18 have different raw ceilings (3.0 and
    # 1.5); the terminator, whose ceiling is between them, rides the Eq. 20
    # weight from one interval to the other -- unless a caller supplied
    # cf_norm_trm, in which case CF_trm is normalized on that instead.
    cf_day = normalize(raw["cf_day_raw"], c.cf_norm_day)
    cf_ngt = normalize(raw["cf_ngt_raw"], c.cf_norm_ngt)
    if c.cf_norm_trm is None:
        cf_trm = normalize_interp(raw["cf_trm_raw"], c.cf_norm_day, c.cf_norm_ngt, b_ngt_trm)
    else:
        cf_trm = normalize(raw["cf_trm_raw"], c.cf_norm_trm)

    cf_comb = blend_confidence(cf_day, cf_trm, cf_ngt, b_ngt_trm, b_trm_day)

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
