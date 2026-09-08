"""CF_cal: the confidence factor on a precision scale instead of a physical one.

CF_comb answers "how much of this branch's ceiling did the dust tests reach"
(Eqs. 16-19, each branch normalized to its own raw ceiling). That is the right
answer for imagery -- 1.0 means saturated whatever the sun is doing, so dust
does not fade at dusk -- but it is not the answer a threshold wants. Night
detection is genuinely harder than day detection, so the same CF_comb bought
worse precision at night: on 42 dust days over East Asia, CF_comb >= 0.2 came
with a false-alarm ratio of 0.399 by day and 0.523 by night.

The two cannot be reconciled by moving the Eq. 19 interval. An interval is two
numbers: it can put the working points where the precision matches, or it can
map the raw ceiling to 1.0, not both (pushing the ceiling up to match precision
caps CF_ngt near 0.55, which is the fading-at-dusk defect again). A monotone
table has no such limit, so this module supplies a second field rather than
bending the first: CF_cal is CF_comb's branches passed through per-branch
monotone tables fitted so that a threshold means the same precision by night
as by day, then blended by Eq. 22 exactly as before. CF_comb is untouched.

What a fitted table can and cannot promise is in :class:`CalibrationTable`.
"""

from dataclasses import dataclass

import numpy as np
import xarray as xr

from shachen.confidence import blend_confidence
from shachen.constants import DEFAULTS, ConfidenceConstants
from shachen.norm import normalize, normalize_cos_zenith

_RAW_VARS = ("cf_day_raw", "cf_trm_raw", "cf_ngt_raw")


@dataclass(frozen=True)
class CalibrationTable:
    """A monotone piecewise-linear map from one branch's raw sum to CF_cal.

    ``knots`` are ``(raw, cf)`` pairs with raw strictly increasing and cf
    non-decreasing, starting at cf 0 and ending at cf 1. Between knots the map
    is linear; below the first knot it is 0, above the last it is 1.

    Fitting one needs ground truth and therefore happens outside this package.
    The construction that gives the knots meaning: for a ladder of levels,
    measure the day branch's false-alarm ratio at that level, then find the
    lowest raw cutoff at which this branch reaches the same false-alarm ratio.
    Two properties follow, and both matter when reading CF_cal:

    * **The promise is about a pool, not an occasion.** The knots equalize
      precision over whatever sample they were fitted on. Precision on any one
      dust event can be far from it -- on the 42-day fit behind
      :data:`EAST_ASIA_AHI`, per-event FAR at CF_cal 0.2 ran from 0.17 to 0.94
      around a pooled 0.40.
    * **Above some level there is nothing to calibrate against.** A branch that
      cannot reach the day branch's precision at high thresholds has no cutoff
      to pin, so the table above its last calibrated knot is monotone only, and
      a threshold there is a ranking, not a precision.
    """

    knots: tuple[tuple[float, float], ...]
    #: Highest cf value that was pinned to a measured precision; above it the
    #: table only ranks. Informational: nothing in the arithmetic reads it.
    calibrated_to: float = 1.0

    def __post_init__(self) -> None:
        if len(self.knots) < 2:
            raise ValueError("a calibration table needs at least two knots")
        raw = np.array([k[0] for k in self.knots], dtype=float)
        cf = np.array([k[1] for k in self.knots], dtype=float)
        if not np.all(np.diff(raw) > 0):
            raise ValueError(f"knot raw values must strictly increase, got {raw.tolist()}")
        if not np.all(np.diff(cf) >= 0):
            raise ValueError(f"knot cf values must not decrease, got {cf.tolist()}")
        if cf[0] != 0.0 or cf[-1] != 1.0:
            raise ValueError(f"a table must run from cf 0 to cf 1, got {cf[0]} to {cf[-1]}")

    @property
    def raw(self) -> np.ndarray:
        return np.array([k[0] for k in self.knots], dtype=float)

    @property
    def cf(self) -> np.ndarray:
        return np.array([k[1] for k in self.knots], dtype=float)

    def __call__(self, raw):
        """Apply the table. NaN in, NaN out; DataArrays keep their coords."""
        values = np.interp(np.asarray(raw, dtype=float), self.raw, self.cf)
        values = np.where(np.isfinite(np.asarray(raw, dtype=float)), values, np.nan)
        if isinstance(raw, xr.DataArray):
            return raw.copy(data=values)
        return values


@dataclass(frozen=True)
class Calibration:
    """The tables CF_cal needs: one per branch that is not the reference.

    The day branch has none because it *is* the reference -- CF_cal equals
    CF_day wherever Eq. 22 is reading the day branch alone.
    """

    night: CalibrationTable
    terminator: CalibrationTable


#: Fitted on 42 dust days over East Asia (Himawari AHI, 19 events, hourly,
#: all hours), against hourly station PM10 >= 600 ug/m3, at a 5x5 pixel median
#: box. Leave-one-event-out, the day-night FAR gap at CF_cal 0.1 / 0.2 / 0.3
#: closes from +0.114 / +0.125 / +0.125 to +0.030 / +0.005 / +0.076, and night
#: AUC rises 0.623 -> 0.633. The night branch stops being calibrated above
#: 0.45: night precision does not reach what the day branch has above that, so
#: there is no cutoff to pin. See the caveats on :class:`CalibrationTable` --
#: especially that the equality is a pooled one.
EAST_ASIA_AHI = Calibration(
    night=CalibrationTable(
        knots=(
            (0.1500, 0.00),
            (0.2985, 0.05),
            (0.3602, 0.10),
            (0.6571, 0.15),
            (0.7411, 0.20),
            (0.7924, 0.25),
            (0.8408, 0.30),
            (0.8877, 0.35),
            (0.9134, 0.40),
            (0.9238, 0.45),
            (1.5000, 1.00),
        ),
        calibrated_to=0.45,
    ),
    terminator=CalibrationTable(
        knots=(
            (0.1500, 0.00),
            (0.2829, 0.05),
            (0.3420, 0.10),
            (1.0694, 0.15),
            (1.3233, 0.20),
            (1.3489, 0.25),
            (1.4543, 0.30),
            (1.5084, 0.35),
            (1.5628, 0.40),
            (2.1088, 0.45),
            (2.1138, 0.50),
            (2.1161, 0.55),
            (2.1322, 0.60),
            (2.1383, 0.65),
            (2.1466, 0.70),
            (2.1478, 0.95),
            (2.5000, 1.00),
        ),
        calibrated_to=0.95,
    ),
)


def calibrate(
    raw: xr.Dataset,
    zenith_deg,
    calibration: Calibration = EAST_ASIA_AHI,
    constants: ConfidenceConstants = DEFAULTS.confidence,
) -> xr.Dataset:
    """CF_cal on the scene grid, from the raw Eqs. 16-18 sums.

    ``raw`` is what :func:`shachen.confidence.confidence_raw` returns. The day
    branch is normalized by Eq. 19 as usual (it is the reference the tables
    were fitted against); the terminator and night branches go through their
    tables; Eq. 22 then blends the three with the same weights CF_comb uses,
    so the two fields agree wherever the sun does not matter.

    Returns ``cf_cal`` plus the three branches behind it (``cf_day``,
    ``cf_trm_cal``, ``cf_ngt_cal``) and the blend weights.
    """
    missing = [name for name in _RAW_VARS if name not in raw]
    if missing:
        raise ValueError(f"raw is missing required variable(s): {missing}")

    b_ngt_trm = normalize_cos_zenith(
        zenith_deg, constants.ngt_trm_zenith_deg, constants.blend_exponent
    )
    b_trm_day = normalize_cos_zenith(
        zenith_deg, constants.trm_day_zenith_deg, constants.blend_exponent
    )
    cf_day = normalize(raw["cf_day_raw"], constants.cf_norm_day)
    cf_trm = calibration.terminator(raw["cf_trm_raw"])
    cf_ngt = calibration.night(raw["cf_ngt_raw"])
    return xr.Dataset(
        {
            "cf_cal": blend_confidence(cf_day, cf_trm, cf_ngt, b_ngt_trm, b_trm_day),
            "cf_day": cf_day,
            "cf_trm_cal": cf_trm,
            "cf_ngt_cal": cf_ngt,
            "b_ngt_trm": b_ngt_trm,
            "b_trm_day": b_trm_day,
        }
    )
