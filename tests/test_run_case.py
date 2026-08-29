"""Acceptance tests for the scripts/run_case.py CLI.

The CLI lives outside the package, so it is loaded by file path. Error paths
run on synthetic/empty directories; the end-to-end test reuses the real
2017-03-23 case (same skipif guard as test_pipeline) and checks the yellow
coloration in the plume box (mean CF 0.48 there gives
red - blue of roughly (1 - D)*CF/1.2, about 0.36, far above the 0.05
threshold; the Gulf control box had mean CF 0.00).
"""

import datetime as dt
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = _ROOT / "scripts" / "run_case.py"

_DATA = _ROOT / "data"
_ABI = _DATA / "2017-03-23-swus" / "abi"
_MERRA = _DATA / "merra2" / "merra2_ts_20170323.nc"
_EMIS = _DATA / "emissivity" / "CAM5K30EM_201703.nc"
_HAVE_CASE = (
    _ABI.exists() and len(list(_ABI.glob("*.nc"))) >= 8 and _MERRA.exists() and _EMIS.exists()
)


def _load_module():
    spec = importlib.util.spec_from_file_location("run_case", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    # Register before exec: the dataclass decorator resolves string
    # annotations through sys.modules[cls.__module__].
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


run_case = _load_module()


# --- input resolution --------------------------------------------------------


def test_unknown_case_raises_keyerror():
    with pytest.raises(KeyError):
        run_case.resolve_inputs("1999-01-01-nowhere")


def test_missing_data_points_at_fetch_case(tmp_path):
    with pytest.raises(FileNotFoundError, match="fetch_case"):
        run_case.resolve_inputs("2017-03-23-swus", data_dir=tmp_path)


def test_cli_rejects_unknown_case_and_color(tmp_path):
    with pytest.raises(SystemExit):
        run_case.main(["1999-01-01-nowhere"])
    with pytest.raises(SystemExit):
        run_case.main(["2017-03-23-swus", "--color", "magenta"])


def test_tuning_registry():
    # The CLI's --tuning presets: "abi" (the default) is the ABI_TUNED
    # retune; "paper" reverts to the published constants.
    from shachen.constants import ABI_TUNED, DEFAULTS

    assert run_case.TUNING["abi"] is ABI_TUNED
    assert run_case.TUNING["paper"] is DEFAULTS


def test_cli_rejects_unknown_tuning():
    with pytest.raises(SystemExit):
        run_case.main(["2017-03-23-swus", "--tuning", "bogus"])


# --- composite background mode ------------------------------------------------


def test_backgrounds_registry():
    assert run_case.BACKGROUNDS == ("semianalytic", "composite")


def test_cli_rejects_unknown_background():
    with pytest.raises(SystemExit):
        run_case.main(["2017-03-23-swus", "--background", "bogus"])


def test_resolve_unknown_background_mode_raises():
    with pytest.raises(ValueError):
        run_case.resolve_inputs("2017-03-23-swus", background="bogus")


def test_resolve_composite_missing_points_at_composite_flag(tmp_path):
    with pytest.raises(FileNotFoundError, match="--composite"):
        run_case.resolve_inputs("2020-12-13-splains", data_dir=tmp_path, background="composite")


# --- output writing (synthetic) ----------------------------------------------


def test_save_outputs_synthetic(tmp_path, latlon_area):
    rng = np.random.default_rng(7)
    ny, nx = latlon_area.shape
    result = xr.Dataset(
        {
            "cf_comb": (("y", "x"), rng.uniform(0.0, 1.0, (ny, nx))),
            "rgb": (
                ("y", "x", "gun"),
                rng.uniform(0.0, 1.0, (ny, nx, 3)),
            ),
            "dust_rgb": (
                ("y", "x", "gun"),
                rng.uniform(0.0, 1.0, (ny, nx, 3)),
            ),
        },
        coords={"gun": ["r", "g", "b"]},
        attrs={
            "area": latlon_area,
            "start_time": dt.datetime(2020, 12, 13, 21, 20),
        },
    )
    nc_path, png_path = run_case.save_outputs(result, tmp_path, "demo")
    assert nc_path.exists() and nc_path.suffix == ".nc"
    assert png_path.exists() and png_path.suffix == ".png"
    assert png_path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    dust_png = tmp_path / "demo_dustrgb.png"
    assert dust_png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    reopened = xr.open_dataset(nc_path)
    assert "cf_comb" in reopened and "rgb" in reopened
    assert isinstance(reopened.attrs["start_time"], str)


# --- end-to-end on the real 2017-03-23 case ----------------------------------


@pytest.mark.skipif(
    not _HAVE_CASE,
    reason="2017-03-23 case incomplete (need ABI + MERRA-2 + CAMEL; "
    "run scripts/fetch_case.py first)",
)
def test_cli_end_to_end_yellow_plume(tmp_path):
    # Explicit --tuning abi (also the default): the plume/control assertions
    # below hold for both constant sets (tuned plume mean 0.452 per the
    # bound sweep documented in docs/deviations.md).
    run_case.main(["2017-03-23-swus", "--outdir", str(tmp_path), "--tuning", "abi"])

    ncs = list(tmp_path.glob("*.nc"))
    pngs = list(tmp_path.glob("*.png"))
    # The DEBRA composite plus the classic Dust RGB baseline PNG.
    assert len(ncs) == 1 and len(pngs) == 2
    assert all(png.stat().st_size > 0 for png in pngs)

    out = xr.open_dataset(ncs[0])
    assert "cf_comb" in out and "rgb" in out
    rgb = np.asarray(out["rgb"].values, dtype=float)
    finite = np.isfinite(rgb)
    assert finite.any()
    assert ((rgb[finite] >= 0.0) & (rgb[finite] <= 1.0)).all()

    from shachen.io.satellite import load_scene

    scene = load_scene(sorted(_ABI.glob("*.nc")))
    lons, lats = scene.attrs["area"].get_lonlats()

    red = rgb[..., 0]
    blu = rgb[..., 2]
    ok = np.isfinite(red) & np.isfinite(blu) & np.isfinite(lons) & np.isfinite(lats)
    plume = ok & (lats >= 30.0) & (lats <= 33.0) & (lons >= -107.0) & (lons <= -104.5)
    control = ok & (lats >= 24.0) & (lats <= 27.0) & (lons >= -90.0) & (lons <= -86.0)
    assert plume.any() and control.any()
    # Yellow dust: blue gun dimmed against red where CF is high, guns nearly
    # equal (grayscale baseline) where CF is near zero.
    assert (red[plume] - blu[plume]).mean() > 0.05
    assert np.abs(red[control] - blu[control]).mean() < 0.02


# --- composite end-to-end on the 2020-12-13 transparent-dust case -------------

_ABI13 = _DATA / "2020-12-13-splains" / "abi"
_COMP13 = _DATA / "2020-12-13-splains" / "composite"
_MERRA13 = _DATA / "merra2" / "merra2_ts_20201213.nc"
_EMIS13 = _DATA / "emissivity" / "CAM5K30EM_202012.nc"


def _composite_day_count() -> int:
    if not _COMP13.exists():
        return 0
    return sum(1 for d in _COMP13.iterdir() if d.is_dir() and len(list(d.glob("*.nc"))) >= 3)


_HAVE_COMPOSITE_CASE = (
    _ABI13.exists()
    and len(list(_ABI13.glob("*.nc"))) >= 8
    and _MERRA13.exists()
    and _EMIS13.exists()
    and _composite_day_count() >= 5
)


@pytest.mark.skipif(
    not _HAVE_COMPOSITE_CASE,
    reason="2020-12-13 composite window not downloaded (run scripts/"
    "fetch_case.py 2020-12-13-splains --composite 14)",
)
def test_composite_recovers_transparent_dust(tmp_path):
    # At the Quick Guide hook the semianalytic background
    # sits above both observed split-window signals, zeroing DT1 and DT2.
    # Measured on the fetched 14-day window: the dusty
    # RSW (-0.99 K) lies below even the observational clear-sky RSW
    # (-0.53 K), so DT1 *cannot* fire at this scan — the structural recovery
    # is DT2 (composite btd_bg -1.57 K vs CAMEL's +0.76 K, obs -1.40 K).
    # The composite run must (1) fire DT2 over a meaningful share of the
    # hook box, (2) raise the box mean confidence over the semianalytic run,
    # and (3) light up a meaningful fraction of the box at all.
    run_case.main(["2020-12-13-splains", "--outdir", str(tmp_path / "semi")])
    run_case.main(
        [
            "2020-12-13-splains",
            "--outdir",
            str(tmp_path / "comp"),
            "--background",
            "composite",
        ]
    )
    (semi_nc,) = (tmp_path / "semi").glob("*.nc")
    (comp_nc,) = (tmp_path / "comp").glob("*.nc")
    assert comp_nc.stem.endswith("_composite")
    semi = xr.open_dataset(semi_nc)
    comp = xr.open_dataset(comp_nc)
    assert "n_valid" in comp

    from shachen.io.satellite import load_scene

    scene = load_scene(sorted(_ABI13.glob("*.nc")))
    lons, lats = scene.attrs["area"].get_lonlats()
    box = (
        np.isfinite(lons)
        & np.isfinite(lats)
        & (lats >= 35.0)
        & (lats <= 36.5)
        & (lons >= -101.5)
        & (lons <= -99.5)
    )
    assert box.any()

    dt2_box = comp["dt2"].values[box]
    dt2_box = dt2_box[np.isfinite(dt2_box)]
    assert dt2_box.size > 0
    assert (dt2_box > 0.05).mean() > 0.10  # measured 0.345

    comp_cf = comp["cf_comb"].values[box]
    comp_cf = comp_cf[np.isfinite(comp_cf)]
    comp_mean = comp_cf.mean()
    semi_mean = np.nanmean(semi["cf_comb"].values[box])
    assert comp_mean > semi_mean  # measured 0.0050 vs 0.0004
    assert (comp_cf > 0.02).mean() > 0.05  # measured 0.101
