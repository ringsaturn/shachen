"""DEBRA dust detection tests DT1-DT3, Eqs. 13-15 (Miller et al. 2017).

DT1/DT2 measure the observed split-window signals against the per-pixel
*dynamic* background (RSW_bg / BTD_bg used as the MIN bound of Eq. 3). DT3 is
the thermal-contrast consensus builder; it is implemented magnitude-reversed
relative to the printed Eq. 15, per the paper's prose ("observations that are
relatively cold compared to MERRA produce high value for DT3"); see
docs/deviations.md:

    DT3 = clip(((T_MERRA - S) - BT10.4) / depth, 0, 1)

with the surface shift S = -10 K (land) / +5 K (ocean) applied as printed
(``T_MERRA - S``) and depth = 50 K.
"""

import numpy as np
import xarray as xr

from shachen.constants import DEFAULTS, DustTestConstants


def _where(condition, if_true, if_false):
    """``np.where`` that preserves xarray labels when any operand is a DataArray."""
    operands = (condition, if_true, if_false)
    if any(isinstance(operand, xr.DataArray) for operand in operands):
        return xr.where(condition, if_true, if_false)
    return np.where(condition, if_true, if_false)


def _dynamic_test(observed, background, maximum: float):
    """Shared Eq. 13/14 kernel: normalize ``observed`` between ``background``.

    The MIN bound of Eq. 3 is the per-pixel dynamic background; MAX is the
    fixed constant. A degenerate range (``background >= maximum``) yields 0.0.
    A NaN background makes the comparison False, so NaN flows through the
    division and propagates as required.
    """
    numerator = observed - background
    denominator = maximum - background
    degenerate = denominator <= 0.0
    # Keep the division well defined where the range is exhausted; those pixels
    # are overwritten with 0.0 below.
    safe_denominator = _where(degenerate, 1.0, denominator)
    value = np.clip(numerator / safe_denominator, 0.0, 1.0)
    return _where(degenerate, 0.0, value)


def dt1(rsw_obs, rsw_bg, constants: DustTestConstants = DEFAULTS.dust_tests):
    """Eq. 13: ``(RSW_obs - RSW_bg) / (MAX_RSW - RSW_bg)`` clipped to [0, 1].

    ``RSW = BT12.3 - BT10.4`` (K). Array-generic. Where the dynamic range is
    degenerate (``RSW_bg >= MAX_RSW``) the test returns 0.0; NaN propagates.
    """
    return _dynamic_test(rsw_obs, rsw_bg, constants.dt1_max_rsw_k)


def dt2(btd_obs, btd_bg, constants: DustTestConstants = DEFAULTS.dust_tests):
    """Eq. 14: like :func:`dt1` for ``BTD = BT8.6 - BT10.4`` with MAX = 3.0 K."""
    return _dynamic_test(btd_obs, btd_bg, constants.dt2_max_btd_k)


def dt3(
    bt_104,
    skin_temperature,
    is_land,
    constants: DustTestConstants = DEFAULTS.dust_tests,
):
    """Eq. 15 (magnitude-reversed, see module docstring).

    ``is_land`` selects the land/ocean surface shift. Array-generic; NaN in
    ``bt_104`` or ``skin_temperature`` propagates.
    """
    shift = _where(is_land, constants.dt3_shift_land_k, constants.dt3_shift_ocean_k)
    reference = skin_temperature - shift
    return np.clip((reference - bt_104) / constants.dt3_depth_k, 0.0, 1.0)


def dust_tests(
    scene: xr.Dataset,
    background: xr.Dataset,
    skin_temperature: xr.DataArray,
    is_land: xr.DataArray,
    constants: DustTestConstants = DEFAULTS.dust_tests,
) -> xr.Dataset:
    """Assemble ``dt1``, ``dt2``, ``dt3`` fields on the scene grid.

    ``scene`` needs ``bt_tir_86``, ``bt_tir_104``, ``bt_tir_123``;
    ``background`` needs ``rsw_bg``, ``btd_bg`` (from
    :func:`shachen.background.background_signals`). All 2-D inputs must share one
    shape (ValueError otherwise). Returns a Dataset with ``dt1``, ``dt2``,
    ``dt3`` in [0, 1] (NaN where inputs are NaN).
    """
    fields = {
        "bt_tir_86": scene["bt_tir_86"],
        "bt_tir_104": scene["bt_tir_104"],
        "bt_tir_123": scene["bt_tir_123"],
        "rsw_bg": background["rsw_bg"],
        "btd_bg": background["btd_bg"],
        "skin_temperature": skin_temperature,
        "is_land": is_land,
    }
    shapes = {name: tuple(field.shape) for name, field in fields.items()}
    if len(set(shapes.values())) > 1:
        raise ValueError(f"dust test inputs must share one shape, got {shapes}")

    bt_104 = fields["bt_tir_104"]
    rsw_obs = fields["bt_tir_123"] - bt_104
    btd_obs = fields["bt_tir_86"] - bt_104
    return xr.Dataset(
        {
            "dt1": dt1(rsw_obs, fields["rsw_bg"], constants),
            "dt2": dt2(btd_obs, fields["btd_bg"], constants),
            "dt3": dt3(bt_104, skin_temperature, is_land, constants),
        }
    )
