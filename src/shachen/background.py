"""Clear-sky background reference signals (Miller et al. 2017, section 3.2).

The semianalytic background: per-band surface emissivity (CAMEL climatology,
interpolated to the sensor band centers) modifies the Planck blackbody radiance
of the MERRA-2 skin temperature, which is inverted back to an equivalent
brightness temperature: ``BT_bg = B^-1(eps * B(T_skin))``. Band differences of
these give the dynamic backgrounds used by the dust tests (Eqs. 13-14):
``RSW_bg = BT_bg(12.3) - BT_bg(10.4)`` and ``BTD_bg = BT_bg(8.6) - BT_bg(10.4)``.

The Planck function is monochromatic at the band centers (``BAND_CENTER_UM``),
with physical constants from scipy.constants; only background *differences*
enter DEBRA, so the monochromatic approximation is well inside the paper's
quoted <5-10% MERRA-driven uncertainty on the final confidence factor.

Missing emissivity (NaN, e.g. ocean in CAMEL) is treated as eps = 1.0: the
paper notes water-surface infrared emissivity is close to unity and needs no
background correction, so the tests fall back to their static bounds there.
"""

import numpy as np
import xarray as xr
from scipy.constants import c as SPEED_OF_LIGHT
from scipy.constants import h as PLANCK_CONSTANT
from scipy.constants import k as BOLTZMANN_CONSTANT

from shachen.constants import BAND_CENTER_UM, Band

#: micrometres -> metres
_UM_TO_M = 1e-6

#: Emissivity used where CAMEL has no retrieval (ocean); paper section 3.2.
_MISSING_EMISSIVITY_FILL = 1.0

#: The two Planck coefficients, as functions of wavelength in metres:
#: c1 = 2 h c^2 (W m2 sr-1), c2 = h c / k (m K).
_C1 = 2.0 * PLANCK_CONSTANT * SPEED_OF_LIGHT**2
_C2 = PLANCK_CONSTANT * SPEED_OF_LIGHT / BOLTZMANN_CONSTANT


def planck_radiance(temperature_k, wavelength_um):
    """Monochromatic Planck spectral radiance B(lambda, T) in W m-2 sr-1 m-1.

    Array-generic (scalars, numpy, xarray); ``wavelength_um`` in micrometres.
    """
    lam = wavelength_um * _UM_TO_M
    return (_C1 / lam**5) / np.expm1(_C2 / (lam * temperature_k))


def planck_temperature(radiance, wavelength_um):
    """Analytic inverse of :func:`planck_radiance`: brightness temperature in K."""
    lam = wavelength_um * _UM_TO_M
    return (_C2 / lam) / np.log1p(_C1 / (lam**5 * radiance))


def _fill_missing_emissivity(emissivity):
    """Replace NaN emissivity with unity, preserving the container type."""
    fillna = getattr(emissivity, "fillna", None)
    if fillna is not None:
        return fillna(_MISSING_EMISSIVITY_FILL)
    return np.where(np.isnan(emissivity), _MISSING_EMISSIVITY_FILL, emissivity)


def background_bt(skin_temperature, emissivity, band: Band):
    """Background brightness temperature ``B^-1(eps * B(T_skin))`` for one band.

    Array-generic; NaN emissivity is treated as 1.0 (identity: BT_bg = T_skin).
    """
    wavelength_um = BAND_CENTER_UM[band]
    eps = _fill_missing_emissivity(emissivity)
    return planck_temperature(eps * planck_radiance(skin_temperature, wavelength_um), wavelength_um)


def background_signals(skin_temperature: xr.DataArray, emissivity: xr.Dataset) -> xr.Dataset:
    """Per-pixel background signals for the dust tests (Eqs. 13-14).

    ``emissivity`` must carry ``emis_tir_86``, ``emis_tir_104``, ``emis_tir_123``
    (the ``load_band_emissivity`` naming), on the same grid as
    ``skin_temperature`` (K). Raises ValueError if the 2-D shapes differ.

    Returns a Dataset with ``rsw_bg`` (= bt_bg_tir_123 - bt_bg_tir_104),
    ``btd_bg`` (= bt_bg_tir_86 - bt_bg_tir_104), and the per-band
    ``bt_bg_tir_86``, ``bt_bg_tir_104``, ``bt_bg_tir_123`` for debugging.
    """
    bands = (Band.TIR_86, Band.TIR_104, Band.TIR_123)

    bt: dict[Band, xr.DataArray] = {}
    for band in bands:
        name = f"emis_{band.value}"
        if name not in emissivity:
            raise ValueError(f"emissivity is missing variable {name!r}")
        eps = emissivity[name]
        if eps.shape != skin_temperature.shape:
            raise ValueError(
                f"shape mismatch: skin_temperature {skin_temperature.shape} vs {name} {eps.shape}"
            )
        bt[band] = background_bt(skin_temperature, eps, band)

    out = xr.Dataset(
        {f"bt_bg_{band.value}": bt[band].rename(f"bt_bg_{band.value}") for band in bands}
    )
    out["rsw_bg"] = bt[Band.TIR_123] - bt[Band.TIR_104]
    out["btd_bg"] = bt[Band.TIR_86] - bt[Band.TIR_104]

    for name, band in (
        ("bt_bg_tir_86", Band.TIR_86),
        ("bt_bg_tir_104", Band.TIR_104),
        ("bt_bg_tir_123", Band.TIR_123),
    ):
        out[name].attrs.update(
            units="K",
            long_name=f"clear-sky background brightness temperature ({band.value})",
        )
    out["rsw_bg"].attrs.update(
        units="K", long_name="background split-window signal BT(12.3) - BT(10.4)"
    )
    out["btd_bg"].attrs.update(
        units="K", long_name="background brightness temperature difference BT(8.6) - BT(10.4)"
    )
    return out
