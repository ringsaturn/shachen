"""shachen: infrared-channel dust algorithms for geostationary imagers.

DEBRA-Dust, the Dynamic Enhancement with Background Reduction Algorithm of
Miller et al. (2017), doi:10.1002/2017JD027365 (with the 26 Feb 2020
erratum), is the primary one; the classic EUMETSAT Dust RGB ships alongside
it as the baseline to compare against.
"""

from shachen.constants import ABI_BANDS, AHI_BANDS, DEBRA_BANDS, DEFAULTS, Band
from shachen.norm import normalize
from shachen.pipeline import run_debra, run_dust_rgb

__all__ = [
    "DEFAULTS",
    "ABI_BANDS",
    "AHI_BANDS",
    "DEBRA_BANDS",
    "Band",
    "normalize",
    "run_debra",
    "run_dust_rgb",
]
