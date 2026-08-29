"""Classic Dust RGB — the naive baseline DEBRA is judged against.

Not part of the DEBRA algorithm (Eqs. 1-29): no background reduction, no cloud
mask, no confidence field, just three fixed stretches of the infrared channels.
RED = BT12.3 - BT10.4, GRN = BT11.2 - BT8.6 (gamma 2.5), BLU = BT10.4. The
stretch values are **per sensor** — the scheme was tuned for SEVIRI and
re-tuned for each later imager, so there is no single canonical set;
:func:`shachen.pipeline.run_dust_rgb` picks one from the scene's reader.
Computed from
:func:`shachen.io.satellite.load_scene` fields so it shares DEBRA's 2-km grid
and the :mod:`shachen.render` path; dust appears pink to magenta over dark
blue-green surfaces.

:func:`shachen.pipeline.run_dust_rgb` is the entry point that runs this over a
scene, the way :func:`shachen.pipeline.run_debra` runs DEBRA.

Where the recipe comes from (the docs page carries these as links):

* Lensky, I. M., and D. Rosenfeld (2008): Clouds-Aerosols-Precipitation
  Satellite Analysis Tool (CAPSAT). Atmos. Chem. Phys., 8, 6739-6753,
  doi:10.5194/acp-8-6739-2008 -- the SEVIRI RGB suite this scheme comes from.
* EUMeTrain, "Compilation of RGB Recipes" -- the SEVIRI Dust RGB as formalised
  for operations; the source of :data:`shachen.constants.DUST_RGB`.
* NOAA/NASA GOES-R "Quick Guide: Dust RGB" (contributor K. Fuell, NASA SPoRT;
  CIRA/RAMMB) -- the ABI band mix used throughout, and the source of
  :data:`shachen.constants.DUST_RGB_ABI`.
* Berndt, E., N. Elmer, L. Schultz, and A. Molthan (2018): A Methodology to
  Determine Recipe Adjustments for Multispectral Composites Derived from
  Next-Generation Advanced Satellite Imagers. J. Atmos. Oceanic Technol., 35,
  643-664, doi:10.1175/JTECH-D-17-0047.1 -- why an ABI recipe needs stretches
  different from SEVIRI's.
"""

import numpy as np
import xarray as xr

from shachen.constants import DUST_RGB, Band, DustRGBConstants
from shachen.norm import normalize

#: The four bands the recipe reads (11.2 um is the non-DEBRA extra).
DUST_RGB_BANDS: tuple[Band, ...] = (Band.TIR_86, Band.TIR_104, Band.TIR_112, Band.TIR_123)

_REQUIRED = tuple(f"bt_{band.value}" for band in DUST_RGB_BANDS)


def dust_rgb(scene: xr.Dataset, constants: DustRGBConstants = DUST_RGB) -> xr.DataArray:
    """The classic Dust RGB composite of ``scene``, floats in [0, 1].

    ``scene`` needs ``bt_tir_86/104/112/123`` (K); a missing variable raises
    ValueError. Returns dims ``(y, x, gun)`` with coordinate
    ``gun = ["r", "g", "b"]`` — the same layout as
    :func:`shachen.imagery.enhanced_rgb`, so :func:`shachen.imagery.to_uint8` and
    :func:`shachen.render.render_debra_png` apply unchanged. NaN propagates.
    """
    missing = [name for name in _REQUIRED if name not in scene]
    if missing:
        raise ValueError(f"scene is missing required variable(s): {missing}")

    bt86 = scene["bt_tir_86"]
    bt104 = scene["bt_tir_104"]
    bt112 = scene["bt_tir_112"]
    bt123 = scene["bt_tir_123"]

    c = constants
    red = normalize(bt123 - bt104, c.red)
    grn = normalize(bt112 - bt86, c.green) ** (1.0 / c.green_gamma)
    blu = normalize(bt104, c.blue)

    return xr.DataArray(
        np.stack([np.asarray(g.values, dtype=float) for g in (red, grn, blu)], axis=-1),
        dims=tuple(bt104.dims) + ("gun",),
        coords={**dict(bt104.coords), "gun": ["r", "g", "b"]},
    )
