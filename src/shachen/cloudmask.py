"""DEBRA cloud mask, Eqs. 1-12 (Miller et al. 2017).

Continuous cloud confidence in [0, 1] from IR-only tests, with dust-restoral
terms so dust pixels are not masked as cloud. Two deliberate deviations from
the *printed* equations, both following the paper's own prose
(the 2020 erratum does not touch either); see docs/deviations.md:

- Eq. 4 (CM2) is implemented magnitude-reversed,
  ``CM2 = 1 - N(BT10.4 - BT6.2; 0, 25)`` (deep convection -> 1, clear -> 0).
  As printed, clear sky (split of 25-55 K) would give CM2 ~ 1 and saturate the
  mask everywhere, degenerating the algorithm; the prose ("in a similar
  fashion" to the reversed Eq. 1) and Figure 3c confirm the reversed intent.
- Eq. 11 (CM_day) uses CM3, not the misprinted CM4 (the 3.9-um test is
  night-only): ``CM_day = (CM1 + CM2 + CM3) * (1 - max(R1, R2_day))``.
"""

import numpy as np
import xarray as xr

from shachen.constants import DEFAULTS, Bounds, CloudMaskConstants
from shachen.norm import normalize

#: Scene variables consumed by the cloud mask (all brightness temperatures, K).
_REQUIRED_BANDS = (
    "bt_swir_39",
    "bt_wv_62",
    "bt_tir_86",
    "bt_tir_104",
    "bt_tir_123",
)


def _check_shapes(scene: xr.Dataset, skin_temperature: xr.DataArray) -> None:
    """ValueError unless every band shares the 2-D shape of ``skin_temperature``."""
    skin_shape = tuple(skin_temperature.shape)
    for name in _REQUIRED_BANDS:
        if name not in scene:
            raise ValueError(f"scene is missing required variable {name!r}")
        band_shape = tuple(scene[name].shape)
        if band_shape != skin_shape:
            raise ValueError(
                f"shape mismatch: scene[{name!r}] has shape {band_shape}, "
                f"skin_temperature has shape {skin_shape}"
            )


def cloud_mask(
    scene: xr.Dataset,
    skin_temperature: xr.DataArray,
    constants: CloudMaskConstants = DEFAULTS.cloud_mask,
) -> xr.Dataset:
    """Cloud-mask components and combined masks on the scene grid.

    ``scene`` needs ``bt_swir_39``, ``bt_wv_62``, ``bt_tir_86``, ``bt_tir_104``,
    ``bt_tir_123`` (K); ``skin_temperature`` (K) must be on the same grid
    (ValueError on 2-D shape mismatch). NaN inputs propagate.

    Returns a Dataset with ``cm1``..``cm4`` (Eqs. 1, 4-6; CM1 reduces to
    ``1 - N(BT10.4 - T_skin; -offset, 0)`` since the per-pixel bounds share a
    constant width), ``r1``, ``r2_day``, ``r2_ngt`` (Eqs. 7-9), ``cm_day``,
    ``cm_ngt`` (Eqs. 10-11), and ``cm_norm_day``, ``cm_norm_ngt`` (Eq. 12).
    """
    c = constants
    _check_shapes(scene, skin_temperature)

    bt39 = scene["bt_swir_39"]
    bt62 = scene["bt_wv_62"]
    bt86 = scene["bt_tir_86"]
    bt104 = scene["bt_tir_104"]
    bt123 = scene["bt_tir_123"]

    # Eqs. 1-2: cold-relative-to-skin test. The per-pixel bounds
    # (T_skin - offset, T_skin) have a constant width, so subtracting T_skin
    # turns them into the fixed bounds (-offset, 0) on the difference.
    cm1_bounds = Bounds(-c.cm1_cold_offset_k, 0.0)
    cm1 = 1.0 - normalize(bt104 - skin_temperature, cm1_bounds)

    # Eq. 4 (magnitude-reversed): deep convection -> 1, clear sky -> 0.
    cm2 = 1.0 - normalize(bt104 - bt62, c.cm2)

    # Eq. 5: day/night split-window thin-cirrus test.
    cm3 = normalize(bt104 - bt123, c.cm3)

    # Eq. 6: night-only 3.9-um thin-cirrus test.
    cm4 = normalize(bt39 - bt104, c.cm4)

    # Eq. 7 (erratum): reverse split window restoral, damped by cloud coldness.
    r1 = normalize(bt123 - bt104, c.r1) * (1.0 - cm1)

    # Eqs. 8-9: 8.6-um restorals, day weighted by CM3, night by CM4.
    r2_base = normalize(bt86 - bt104, c.r2)
    r2_day = r2_base * (1.0 - cm2) * (1.0 - cm3)
    r2_ngt = r2_base * (1.0 - cm2) * (1.0 - cm4)

    # Eqs. 10-11 (Eq. 11 uses CM3; the printed CM4 is a typo).
    cm_ngt = (cm1 + cm2 + cm4) * (1.0 - np.maximum(r1, r2_ngt))
    cm_day = (cm1 + cm2 + cm3) * (1.0 - np.maximum(r1, r2_day))

    # Eq. 12.
    cm_norm_day = normalize(cm_day, c.cm_norm)
    cm_norm_ngt = normalize(cm_ngt, c.cm_norm)

    return xr.Dataset(
        {
            "cm1": cm1,
            "cm2": cm2,
            "cm3": cm3,
            "cm4": cm4,
            "r1": r1,
            "r2_day": r2_day,
            "r2_ngt": r2_ngt,
            "cm_day": cm_day,
            "cm_ngt": cm_ngt,
            "cm_norm_day": cm_norm_day,
            "cm_norm_ngt": cm_norm_ngt,
        }
    )
