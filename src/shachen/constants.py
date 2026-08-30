"""All frozen numbers in one place: DEBRA's tuning constants and the
classic Dust RGB recipe.

Every (MIN, MAX) scaling bound, offset, and weight from Miller et al. (2017),
JGR Atmospheres, doi:10.1002/2017JD027365, using the 26 Feb 2020 erratum
versions of Eqs. 7, 21, 22, 24, 25. The paper notes "minor retuning" may be
needed per sensor; retune here, nowhere else.

Sign/orientation conventions follow the printed equations exactly.

The sensor facts here (:class:`Band` and the channel maps) are not from the
paper; they are what each instrument calls the roles the algorithms read.
"""

from dataclasses import dataclass, field
from enum import StrEnum


@dataclass(frozen=True)
class Bounds:
    """(MIN, MAX) pair for the Eq. 3 normalization primitive."""

    min: float
    max: float


class Band(StrEnum):
    """Spectral roles the algorithms read (nominal wavelengths in um).

    Seven of these are DEBRA's inputs (:data:`DEBRA_BANDS`);
    ``TIR_112`` exists only for the classic Dust RGB green gun.
    """

    VIS_064 = "vis_064"  #: cloud mask + day baseline image
    NIR_160 = "nir_160"  #: daytime cloud test (reserved)
    SWIR_39 = "swir_39"  #: night thin-cirrus test CM4
    WV_62 = "wv_62"  #: deep-convection test CM2
    TIR_86 = "tir_86"  #: dust test DT2 (8.4-8.6 um)
    TIR_104 = "tir_104"  #: clean window reference (10.3-10.4 um)
    TIR_112 = "tir_112"  #: classic Dust RGB green gun only (11.2 um; not a DEBRA input)
    TIR_123 = "tir_123"  #: dirty window, RSW / DT1 (12.3 um)


#: Spectral role -> GOES-R ABI channel name (satpy dataset names).
ABI_BANDS: dict[Band, str] = {
    Band.VIS_064: "C02",
    Band.NIR_160: "C05",
    Band.SWIR_39: "C07",
    Band.WV_62: "C08",
    Band.TIR_86: "C11",
    Band.TIR_104: "C13",
    Band.TIR_112: "C14",
    Band.TIR_123: "C15",
}

#: Spectral role -> Himawari AHI band name (satpy dataset names).
AHI_BANDS: dict[Band, str] = {
    Band.VIS_064: "B03",
    Band.NIR_160: "B05",
    Band.SWIR_39: "B07",
    Band.WV_62: "B08",
    Band.TIR_86: "B11",
    Band.TIR_104: "B13",
    Band.TIR_112: "B14",
    Band.TIR_123: "B15",
}

#: The seven DEBRA algorithm input roles (Miller et al. 2017 Table 1).
#: ``Band.TIR_112`` exists only for the classic Dust RGB comparison baseline
#: (:mod:`shachen.dustrgb`); pipelines that feed DEBRA alone load these.
DEBRA_BANDS: tuple[Band, ...] = tuple(band for band in Band if band is not Band.TIR_112)

#: Nominal central wavelengths (um) used to interpolate the UWBF emissivity
#: hinge points to sensor band centers (ABI values; AHI within tolerance).
BAND_CENTER_UM: dict[Band, float] = {
    Band.VIS_064: 0.64,
    Band.NIR_160: 1.61,
    Band.SWIR_39: 3.90,
    Band.WV_62: 6.19,
    Band.TIR_86: 8.44,
    Band.TIR_104: 10.33,
    Band.TIR_112: 11.19,
    Band.TIR_123: 12.30,
}


@dataclass(frozen=True)
class CloudMaskConstants:
    """Eqs. 1-12."""

    #: CM1 (Eq. 1): 1 - N(BT10.4; T_skin - cm1_cold_offset_k, T_skin)
    cm1_cold_offset_k: float = 50.0
    #: CM2 (Eq. 4): N(BT10.4 - BT6.2)
    cm2: Bounds = field(default_factory=lambda: Bounds(0.0, 25.0))
    #: CM3 (Eq. 5): N(BT10.4 - BT12.3), day/night thin cirrus
    cm3: Bounds = field(default_factory=lambda: Bounds(2.0, 4.5))
    #: CM4 (Eq. 6): N(BT3.9 - BT10.4), night-only thin cirrus
    cm4: Bounds = field(default_factory=lambda: Bounds(5.0, 8.0))
    #: R1 (Eq. 7, erratum): N(BT12.3 - BT10.4) * (1 - CM1)
    r1: Bounds = field(default_factory=lambda: Bounds(0.0, 3.5))
    #: R2 (Eqs. 8-9): N(BT8.6 - BT10.4) * restoral weights
    r2: Bounds = field(default_factory=lambda: Bounds(-1.0, 3.0))
    #: Eq. 12: CM_norm = N(CM)
    cm_norm: Bounds = field(default_factory=lambda: Bounds(0.45, 0.80))


@dataclass(frozen=True)
class DustTestConstants:
    """Eqs. 13-15."""

    #: DT1 (Eq. 13): (RSW_obs - RSW_bg) / (max_rsw_k - RSW_bg), RSW = BT12.3 - BT10.4
    dt1_max_rsw_k: float = 3.5
    #: DT2 (Eq. 14): (BTD_obs - BTD_bg) / (max_btd_k - BTD_bg), BTD = BT8.6 - BT10.4
    dt2_max_btd_k: float = 3.0
    #: DT3 (Eq. 15): (BT10.4 - (T_merra - S - dt3_depth_k)) / dt3_depth_k
    dt3_shift_land_k: float = -10.0
    dt3_shift_ocean_k: float = 5.0
    dt3_depth_k: float = 50.0


@dataclass(frozen=True)
class ConfidenceConstants:
    """Eqs. 16-22 (Eqs. 21-22 per erratum)."""

    #: Eq. 17: CF_trm weights DT3 by this factor; Eq. 18: CF_ngt likewise.
    dt3_weight_trm: float = 0.5
    dt3_weight_ngt: float = 0.5
    #: Eq. 19: each CF variant normalized with these bounds.
    cf_norm: Bounds = field(default_factory=lambda: Bounds(0.25, 2.50))
    #: Eqs. 20-21: terminator blending on cos(theta_sun), exponent 1.5.
    blend_exponent: float = 1.5
    #: night/terminator interface: N(cos theta; cos 105 deg, cos 90 deg)
    ngt_trm_zenith_deg: Bounds = field(default_factory=lambda: Bounds(105.0, 90.0))
    #: terminator/day interface: N(cos theta; cos 90 deg, cos 75 deg)
    trm_day_zenith_deg: Bounds = field(default_factory=lambda: Bounds(90.0, 75.0))


@dataclass(frozen=True)
class ImageryConstants:
    """Eqs. 23-29 (Eqs. 24-25 per erratum)."""

    #: Eq. 25 (erratum): B_bg = 1 - N(theta_sun; 79 deg, 89 deg)^1.5
    bg_blend_zenith_deg: Bounds = field(default_factory=lambda: Bounds(79.0, 89.0))
    bg_blend_exponent: float = 1.5
    #: Eqs. 27-29: RED = GRN = BI*(1 - min(CF, cf_cap)) + CF;
    #:             BLU = BI*(1 - min(CF, cf_cap)) + blue_dimming * CF
    cf_cap: float = 0.5
    blue_dimming: float = 0.10
    #: Per-gun rescale [0, gun_max] -> [0, 255]
    gun_max: float = 1.2


#: Per-gun dimming triples (D_red, D_grn, D_blu) generalizing Eqs. 27-29: each
#: gun is BI*(1 - min(CF, cf_cap)) + D_gun*CF, a full gun has D = 1.0. Presets
#: from the paper's section 4.2 text ("dust can be made pink by applying
#: D = 0.25 to equations (28) and (29), green by applying D = 0.10 to
#: equations (27) and (29), or blue by applying D = 0.25 to (27) and (28)").
COLOR_DIMMING: dict[str, tuple[float, float, float]] = {
    "yellow": (1.0, 1.0, 0.10),
    "pink": (1.0, 0.25, 0.25),
    "green": (0.10, 1.0, 0.10),
    "blue": (0.25, 0.25, 1.0),
}


#: Cloud-cleared composite background (paper section 3.2, Option B): rolling
#: window length ("multiday (e.g., 14 day) composites") and the fewest
#: composite days accepted at load time before the background is refused.
COMPOSITE_WINDOW_DAYS: int = 14
COMPOSITE_MIN_DAYS: int = 5


@dataclass(frozen=True)
class DustRGBConstants:
    """One classic Dust RGB stretch set. Not part of DEBRA (Eqs. 1-29).

    Each gun is ``N(x; MIN, MAX) ** (1/gamma)`` over the GOES-R/Himawari band
    mix of the CIRA Dust RGB Quick Guide (satpy's ``dust`` composite):
    12.3-10.3, 11.2-8.4, 10.3 — note the green minuend is 11.2 um, not the
    10.8 um of SEVIRI's ``IR10.8 - IR8.7``.

    The stretch values are **per sensor**: :data:`DUST_RGB` (SEVIRI) and
    :data:`DUST_RGB_ABI` are the two published sets, dispatched by
    :data:`DUST_RGB_BY_READER`. See :mod:`shachen.dustrgb` for the references.
    """

    #: RED: BT12.3 - BT10.4 (split window), [-4, +2] K, gamma 1.
    red: Bounds = field(default_factory=lambda: Bounds(-4.0, 2.0))
    #: GRN: BT11.2 - BT8.6, [0, 15] K, gamma 2.5 (note: 11.2 um, not 10.4).
    green: Bounds = field(default_factory=lambda: Bounds(0.0, 15.0))
    green_gamma: float = 2.5
    #: BLU: BT10.4, [261, 289] K, gamma 1.
    blue: Bounds = field(default_factory=lambda: Bounds(261.0, 289.0))


#: The original SEVIRI stretches (Lensky and Rosenfeld 2008, as formalised in
#: EUMeTrain's recipe compilation; satpy's generic ``dust_default``
#: enhancement). Applied to Himawari here, which is what satpy does too: it
#: ships no ``dust_ahi`` override. Also the fallback for an unknown sensor.
DUST_RGB = DustRGBConstants()

#: The ABI-adjusted stretches from the CIRA GOES-R Dust RGB Quick Guide (satpy's
#: sensor-specific ``dust_abi`` enhancement), converted from the Quick Guide's
#: degrees Celsius. Berndt et al. (2018) is the method behind the adjustment:
#: ABI's channels do not sit where SEVIRI's do, so the SEVIRI numbers put the
#: clear-sky land background at the wrong end of the red stretch.
DUST_RGB_ABI = DustRGBConstants(
    red=Bounds(-6.7, 2.6),
    green=Bounds(-0.5, 20.0),
    blue=Bounds(261.2, 288.7),
)

#: satpy reader -> the stretch set that sensor's operational product uses.
#: :func:`shachen.pipeline.run_dust_rgb` reads this off ``scene.attrs["reader"]``
#: so a baseline rendered here matches the image forecasters see. The same
#: rule is why there is no single canonical recipe to freeze: the scheme is
#: retuned per imager (satpy carries ``ash_abi``, ``convection_abi`` and the
#: night-microphysics variants for the same reason).
DUST_RGB_BY_READER: dict[str, DustRGBConstants] = {
    "abi_l1b": DUST_RGB_ABI,
    "ahi_hsd": DUST_RGB,
}


@dataclass(frozen=True)
class DebraConstants:
    """The whole tuning surface, one group per stage of the algorithm."""

    cloud_mask: CloudMaskConstants = field(default_factory=CloudMaskConstants)
    dust_tests: DustTestConstants = field(default_factory=DustTestConstants)
    confidence: ConfidenceConstants = field(default_factory=ConfidenceConstants)
    imagery: ImageryConstants = field(default_factory=ImageryConstants)


#: The paper-tuned default constant set.
DEFAULTS = DebraConstants()

#: Optional retune for the ABI + MERRA-2 + CAMEL stack. Single deviation from
#: the paper: the Eq. 19 lower bound is raised 0.25 -> 0.40. Rationale: on this
#: ancillary stack DT3 carries a 0.2-0.5 clear-sky floor over land
#: (MERRA-2 skin T runs warmer than BT10.4, plus the -10 K land shift), which
#: leaks through low-cloud-mask pixels as a faint yellow tint over vegetated
#: and low-cloud areas. Sweeping the bound over the three reference cases:
#: 0.40 removes 92% of the tinted area (SE-US box, fraction CF > 0.05:
#: 0.37 -> 0.028) while the 2017-03-23 plume mean CF drops only 5%
#: (0.476 -> 0.452) and the 2020-12-23 plume core stays distinct (p90 0.40).
#: The paper itself anticipates "minor retuning" per sensor.
ABI_TUNED = DebraConstants(
    confidence=ConfidenceConstants(cf_norm=Bounds(0.40, 2.50)),
)
