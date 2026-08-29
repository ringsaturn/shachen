"""Smoke test for the L1b scene loader.

Runs only when the 2017-03-23 case has been fetched
(``scripts/fetch_case.py 2017-03-23-swus --no-ancillary``).
"""

from pathlib import Path

import pytest

CASE_DIR = Path(__file__).resolve().parents[1] / "data" / "2017-03-23-swus" / "abi"

pytestmark = pytest.mark.skipif(
    len(list(CASE_DIR.glob("*.nc"))) < 8 if CASE_DIR.exists() else True,
    reason="2017-03-23 case not downloaded or missing C14 (re-run scripts/fetch_case.py)",
)


@pytest.fixture(scope="module")
def scene():
    from shachen.io.satellite import load_scene

    return load_scene(sorted(CASE_DIR.glob("*.nc")))


def test_all_bands_on_common_grid(scene):
    assert sorted(scene.data_vars) == [
        "bt_swir_39",
        "bt_tir_104",
        "bt_tir_112",
        "bt_tir_123",
        "bt_tir_86",
        "bt_wv_62",
        "refl_nir_160",
        "refl_vis_064",
    ]
    shapes = {v.shape for v in scene.data_vars.values()}
    assert len(shapes) == 1  # everything on the 2-km grid
    assert scene.attrs["area"].pixel_size_x == pytest.approx(2004, abs=5)


def test_bt_ranges_physical(scene):
    bt = scene["bt_tir_104"]
    assert float(bt.min()) > 180.0
    assert float(bt.max()) < 340.0


def test_split_window_difference_small(scene):
    # BT10.4 - BT12.3 over the scene should be within a few kelvin.
    d = (scene["bt_tir_104"] - scene["bt_tir_123"]).mean()
    assert -5.0 < float(d) < 5.0


def test_tir_only_roles_subset(scene):
    # Composite-background days hold only the three TIR band files
    # and are loaded with roles=COMPOSITE_BANDS.
    from shachen.composite import COMPOSITE_BANDS
    from shachen.io.satellite import load_scene

    tir_files = [
        f
        for f in sorted(CASE_DIR.glob("*.nc"))
        if any(f"{chan}_" in f.name for chan in ("C11", "C13", "C15"))
    ]
    assert len(tir_files) == 3

    tir = load_scene(tir_files, roles=COMPOSITE_BANDS)
    assert sorted(tir.data_vars) == ["bt_tir_104", "bt_tir_123", "bt_tir_86"]
    assert tir["bt_tir_104"].shape == scene["bt_tir_104"].shape
    assert tir.attrs["area"].pixel_size_x == pytest.approx(2004, abs=5)
    assert tir.attrs["start_time"] == scene.attrs["start_time"]
