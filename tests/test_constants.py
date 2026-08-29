"""Check the constants table against values transcribed from the paper.

Each expected value here was read directly from Miller et al. (2017) /
the 2020 erratum, independently of constants.py, so a typo in either
place fails the test.
"""

import pytest

from shachen.constants import (
    ABI_BANDS,
    ABI_TUNED,
    AHI_BANDS,
    BAND_CENTER_UM,
    COLOR_DIMMING,
    DEFAULTS,
    Band,
)


class TestCloudMaskBounds:
    def test_cm1_cold_offset(self):
        # Eq. 1: T_cold = T_skin - 50 K
        assert DEFAULTS.cloud_mask.cm1_cold_offset_k == 50.0

    def test_cm2_deep_convection(self):
        # Eq. 4: N(BT10.4 - BT6.2; 0, 25)
        assert (DEFAULTS.cloud_mask.cm2.min, DEFAULTS.cloud_mask.cm2.max) == (0.0, 25.0)

    def test_cm3_split_window_cirrus(self):
        # Eq. 5: N(BT10.4 - BT12; 2.0, 4.5)
        assert (DEFAULTS.cloud_mask.cm3.min, DEFAULTS.cloud_mask.cm3.max) == (2.0, 4.5)

    def test_cm4_night_cirrus(self):
        # Eq. 6: N(BT3.9 - BT10.4; 5.0, 8.0)
        assert (DEFAULTS.cloud_mask.cm4.min, DEFAULTS.cloud_mask.cm4.max) == (5.0, 8.0)

    def test_r1_dust_restoral(self):
        # Eq. 7 (erratum): N(BT12 - BT10.4; 0.0, 3.5)
        assert (DEFAULTS.cloud_mask.r1.min, DEFAULTS.cloud_mask.r1.max) == (0.0, 3.5)

    def test_r2_dust_restoral(self):
        # Eqs. 8-9: N(BT8.6 - BT10.4; -1.0, 3.0)
        assert (DEFAULTS.cloud_mask.r2.min, DEFAULTS.cloud_mask.r2.max) == (-1.0, 3.0)

    def test_cm_norm(self):
        # Eq. 12: N(CM; 0.45, 0.80)
        assert (DEFAULTS.cloud_mask.cm_norm.min, DEFAULTS.cloud_mask.cm_norm.max) == (
            0.45,
            0.80,
        )


class TestDustTestBounds:
    def test_dt1_max_rsw(self):
        # Eq. 13: MAX_RSW = 3.5 K
        assert DEFAULTS.dust_tests.dt1_max_rsw_k == 3.5

    def test_dt2_max_btd(self):
        # Eq. 14: MAX(8.6-10.4) = 3.0 K
        assert DEFAULTS.dust_tests.dt2_max_btd_k == 3.0

    def test_dt3_shifts_and_depth(self):
        # Eq. 15: S = -10 K (land), +5 K (ocean); depth 50 K
        assert DEFAULTS.dust_tests.dt3_shift_land_k == -10.0
        assert DEFAULTS.dust_tests.dt3_shift_ocean_k == 5.0
        assert DEFAULTS.dust_tests.dt3_depth_k == 50.0


class TestConfidenceBounds:
    def test_dt3_weights(self):
        # Eqs. 17-18: 0.5 * DT3 in terminator and night CFs
        assert DEFAULTS.confidence.dt3_weight_trm == 0.5
        assert DEFAULTS.confidence.dt3_weight_ngt == 0.5

    def test_cf_norm(self):
        # Eq. 19: N(CF; 0.25, 2.50)
        assert (DEFAULTS.confidence.cf_norm.min, DEFAULTS.confidence.cf_norm.max) == (
            0.25,
            2.50,
        )

    def test_terminator_blend(self):
        # Eqs. 20-21: exponent 1.5; night/trm (105 deg, 90 deg); trm/day (90 deg, 75 deg)
        c = DEFAULTS.confidence
        assert c.blend_exponent == 1.5
        assert (c.ngt_trm_zenith_deg.min, c.ngt_trm_zenith_deg.max) == (105.0, 90.0)
        assert (c.trm_day_zenith_deg.min, c.trm_day_zenith_deg.max) == (90.0, 75.0)


class TestImageryBounds:
    def test_bg_blend(self):
        # Eq. 25 (erratum): B_bg = 1 - N(theta; 79 deg, 89 deg)^1.5
        i = DEFAULTS.imagery
        assert (i.bg_blend_zenith_deg.min, i.bg_blend_zenith_deg.max) == (79.0, 89.0)
        assert i.bg_blend_exponent == 1.5

    def test_color_modulation(self):
        # Eqs. 27-29: min(CF, 0.5); D = 0.10; rescale [0, 1.2] -> [0, 255]
        i = DEFAULTS.imagery
        assert i.cf_cap == 0.5
        assert i.blue_dimming == 0.10
        assert i.gun_max == 1.2


class TestAbiTuned:
    def test_single_deviation_is_cf_norm_floor(self):
        # ABI retune: only the Eq. 19 lower bound moves, 0.25 -> 0.40.
        t = ABI_TUNED.confidence.cf_norm
        assert (t.min, t.max) == (0.40, 2.50)

    def test_everything_else_matches_paper(self):
        assert ABI_TUNED.cloud_mask == DEFAULTS.cloud_mask
        assert ABI_TUNED.dust_tests == DEFAULTS.dust_tests
        assert ABI_TUNED.imagery == DEFAULTS.imagery
        c, d = ABI_TUNED.confidence, DEFAULTS.confidence
        assert c.dt3_weight_trm == d.dt3_weight_trm
        assert c.dt3_weight_ngt == d.dt3_weight_ngt
        assert c.blend_exponent == d.blend_exponent
        assert c.ngt_trm_zenith_deg == d.ngt_trm_zenith_deg
        assert c.trm_day_zenith_deg == d.trm_day_zenith_deg


class TestColorDimming:
    def test_presets_match_paper(self):
        # Section 4.2: yellow D=0.10 on blue; pink D=0.25 on grn+blu;
        # green D=0.10 on red+blu; blue D=0.25 on red+grn.
        assert COLOR_DIMMING["yellow"] == (1.0, 1.0, 0.10)
        assert COLOR_DIMMING["pink"] == (1.0, 0.25, 0.25)
        assert COLOR_DIMMING["green"] == (0.10, 1.0, 0.10)
        assert COLOR_DIMMING["blue"] == (0.25, 0.25, 1.0)

    def test_yellow_matches_blue_dimming(self):
        assert COLOR_DIMMING["yellow"][2] == DEFAULTS.imagery.blue_dimming


class TestBandMappings:
    @pytest.mark.parametrize(
        "band,abi,ahi",
        [
            (Band.VIS_064, "C02", "B03"),
            (Band.NIR_160, "C05", "B05"),
            (Band.SWIR_39, "C07", "B07"),
            (Band.WV_62, "C08", "B08"),
            (Band.TIR_86, "C11", "B11"),
            (Band.TIR_104, "C13", "B13"),
            (Band.TIR_123, "C15", "B15"),
        ],
    )
    def test_mapping(self, band, abi, ahi):
        assert ABI_BANDS[band] == abi
        assert AHI_BANDS[band] == ahi

    def test_all_seven_roles_mapped(self):
        assert set(ABI_BANDS) == set(Band)
        assert set(AHI_BANDS) == set(Band)
        assert set(BAND_CENTER_UM) == set(Band)

    def test_band_centers_ordered(self):
        centers = [BAND_CENTER_UM[b] for b in Band]
        assert centers == sorted(centers)
