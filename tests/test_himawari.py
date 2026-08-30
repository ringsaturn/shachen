"""Acceptance tests for Himawari AHI support.

Network-free tests cover the fetch_case.py AHI additions (case registry,
sensor dispatch tables, timeline prefix, R20 key selection), the sensor-aware
run_case.resolve_inputs, and the load_scene bbox signature. The real cropped
scene load and the end-to-end CLI run are gated on the fetched
2016-04-21-mongolia case.
"""

import dataclasses
import datetime as dt
import importlib.util
import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS_DIR = _ROOT / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import fetch_case  # noqa: E402

from shachen.constants import AHI_BANDS, Band  # noqa: E402

_CASE = "2016-04-21-mongolia"
_WHEN = dt.datetime(2016, 4, 21, 8, 0)
_BBOX = (80.0, 30.0, 135.0, 55.0)

_AHI_NAMES = tuple(
    f"HS_H08_20160421_0800_{chan}_FLDK_R20_S0101.DAT.bz2" for chan in fetch_case.AHI_CHANNELS
)

#: On-disk names: B03 is stored decompressed (satpy's bz2 patterns hardcode
#: native resolutions), the rest as downloaded.
_LOCAL_NAMES = tuple(name[: -len(".bz2")] if "_B03_" in name else name for name in _AHI_NAMES)


def _load_run_case():
    if "run_case" in sys.modules:
        return sys.modules["run_case"]
    spec = importlib.util.spec_from_file_location("run_case", _SCRIPTS_DIR / "run_case.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


run_case = _load_run_case()


# --- case registry and dispatch tables ----------------------------------------


def test_mongolia_case_registered():
    case = fetch_case.CASES[_CASE]
    assert case.bucket == "noaa-himawari8"
    assert case.product == "AHI-L1b-FLDK"
    assert case.when == _WHEN
    assert case.sensor == "ahi"
    assert case.bbox == _BBOX


def test_abi_cases_keep_defaults():
    case = fetch_case.CASES["2017-03-23-swus"]
    assert case.sensor == "abi"
    assert case.bbox is None


def test_ahi_channels_match_band_roles():
    assert fetch_case.AHI_CHANNELS == ("B03", "B05", "B07", "B08", "B11", "B13", "B14", "B15")
    assert fetch_case.AHI_CHANNELS == tuple(AHI_BANDS[band] for band in Band)


def test_sensor_dispatch_tables():
    assert fetch_case.SENSOR_READERS == {"abi": "abi_l1b", "ahi": "ahi_hsd"}
    assert fetch_case.SENSOR_GLOBS == {"abi": "*.nc", "ahi": "*.DAT*"}
    assert fetch_case.SENSOR_CHANNELS == {
        "abi": fetch_case.ABI_CHANNELS,
        "ahi": fetch_case.AHI_CHANNELS,
    }


# --- timeline prefix ----------------------------------------------------------


def test_ahi_timeline_prefix():
    case = fetch_case.CASES[_CASE]
    assert fetch_case.ahi_timeline_prefix(case) == "noaa-himawari8/AHI-L1b-FLDK/2016/04/21/0800/"


def test_ahi_timeline_prefix_rejects_off_timeline_times():
    case = fetch_case.CASES[_CASE]
    with pytest.raises(ValueError):
        fetch_case.ahi_timeline_prefix(dataclasses.replace(case, when=_WHEN.replace(minute=5)))
    with pytest.raises(ValueError):
        fetch_case.ahi_timeline_prefix(dataclasses.replace(case, when=_WHEN.replace(second=30)))


# --- R20 key selection --------------------------------------------------------

_PREFIX = "noaa-himawari8/AHI-L1b-FLDK/2016/04/21/0800/"


def _listing(names: tuple[str, ...]) -> list[str]:
    return [_PREFIX + name for name in names]


def test_select_ahi_keys_picks_r20_per_channel():
    # Full timeline flavor: R05/R10/R20 triplicate for B03, R10/R20 for the
    # 1-km distractor bands, R20 for the 2-km bands.
    names = (
        "HS_H08_20160421_0800_B01_FLDK_R10_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B01_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B02_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B03_FLDK_R05_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B03_FLDK_R10_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B03_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B04_FLDK_R10_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B05_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B07_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B08_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B11_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B13_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B14_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B15_FLDK_R20_S0101.DAT.bz2",
        "HS_H08_20160421_0800_B16_FLDK_R20_S0101.DAT.bz2",
    )
    selected = fetch_case.select_ahi_keys(_listing(names))
    assert selected == _listing(_AHI_NAMES)


def test_select_ahi_keys_multiple_segments_sorted():
    names = (
        "HS_H08_20160421_0800_B03_FLDK_R20_S0210.DAT.bz2",
        "HS_H08_20160421_0800_B03_FLDK_R20_S0110.DAT.bz2",
    )
    selected = fetch_case.select_ahi_keys(_listing(names), channels=("B03",))
    assert selected == [
        _PREFIX + "HS_H08_20160421_0800_B03_FLDK_R20_S0110.DAT.bz2",
        _PREFIX + "HS_H08_20160421_0800_B03_FLDK_R20_S0210.DAT.bz2",
    ]


def test_select_ahi_keys_native_layout_falls_back_to_coarsest():
    # Post-2019 timeline flavor: no R20 repackaging for the reflective
    # bands — B03 native R05 / B05 native R10, multi-segment, unordered.
    names = (
        "HS_H08_20210315_0000_B03_FLDK_R05_S0210.DAT.bz2",
        "HS_H08_20210315_0000_B03_FLDK_R05_S0110.DAT.bz2",
        "HS_H08_20210315_0000_B05_FLDK_R10_S0110.DAT.bz2",
        "HS_H08_20210315_0000_B13_FLDK_R20_S0110.DAT.bz2",
    )
    selected = fetch_case.select_ahi_keys(_listing(names), channels=("B03", "B05", "B13"))
    assert selected == _listing(
        (
            "HS_H08_20210315_0000_B03_FLDK_R05_S0110.DAT.bz2",
            "HS_H08_20210315_0000_B03_FLDK_R05_S0210.DAT.bz2",
            "HS_H08_20210315_0000_B05_FLDK_R10_S0110.DAT.bz2",
            "HS_H08_20210315_0000_B13_FLDK_R20_S0110.DAT.bz2",
        )
    )


def test_select_ahi_keys_missing_channel_raises():
    names = tuple(n for n in _AHI_NAMES if "_B13_" not in n)
    with pytest.raises(RuntimeError, match="B13"):
        fetch_case.select_ahi_keys(_listing(names))


def test_ahi_local_name_decompresses_repackaged_r20_b03_only():
    b03_r20 = "HS_H08_20160421_0800_B03_FLDK_R20_S0101.DAT.bz2"
    b03_r05 = "HS_H08_20210315_0000_B03_FLDK_R05_S0110.DAT.bz2"
    b13 = "HS_H08_20160421_0800_B13_FLDK_R20_S0101.DAT.bz2"
    assert fetch_case.ahi_local_name(_PREFIX + b03_r20) == b03_r20[: -len(".bz2")]
    # Native-resolution B03 matches satpy's compressed pattern: keep the bz2.
    assert fetch_case.ahi_local_name(_PREFIX + b03_r05) == b03_r05
    assert fetch_case.ahi_local_name(_PREFIX + b13) == b13


# --- fetch dispatch -----------------------------------------------------------


def test_fetch_l1b_dispatches_on_sensor(monkeypatch, tmp_path):
    calls: list[str] = []
    monkeypatch.setattr(fetch_case, "fetch_abi", lambda case, out_dir: calls.append("abi") or [])
    monkeypatch.setattr(fetch_case, "fetch_ahi", lambda case, out_dir: calls.append("ahi") or [])
    fetch_case.fetch_l1b(fetch_case.CASES["2017-03-23-swus"], tmp_path)
    fetch_case.fetch_l1b(fetch_case.CASES[_CASE], tmp_path)
    assert calls == ["abi", "ahi"]


def test_fetch_l1b_unknown_sensor_raises(tmp_path):
    bogus = dataclasses.replace(fetch_case.CASES[_CASE], sensor="modis")
    with pytest.raises(ValueError):
        fetch_case.fetch_l1b(bogus, tmp_path)


# --- run_case input resolution ------------------------------------------------


def _make_tree(root: Path, skip_channel: str | None = None) -> None:
    ahi_dir = root / _CASE / "ahi"
    ahi_dir.mkdir(parents=True)
    for name in _LOCAL_NAMES:
        if skip_channel is not None and f"_{skip_channel}_" in name:
            continue
        (ahi_dir / name).touch()
    (root / "merra2").mkdir()
    (root / "merra2" / "merra2_ts_20160421.nc").touch()
    (root / "emissivity").mkdir()
    (root / "emissivity" / "CAM5K30EM_201604.nc").touch()


def test_caseinputs_fields():
    names = {f.name for f in dataclasses.fields(run_case.CaseInputs)}
    assert {"l1b_files", "sensor", "bbox"} <= names
    assert "abi_files" not in names


def test_resolve_mongolia_missing_points_at_fetch_case(tmp_path):
    with pytest.raises(FileNotFoundError, match="fetch_case"):
        run_case.resolve_inputs(_CASE, data_dir=tmp_path)


def test_resolve_mongolia_synthetic_tree(tmp_path):
    _make_tree(tmp_path)
    inputs = run_case.resolve_inputs(_CASE, data_dir=tmp_path)
    assert inputs.sensor == "ahi"
    assert inputs.bbox == _BBOX
    assert inputs.when == _WHEN
    assert len(inputs.l1b_files) == len(fetch_case.AHI_CHANNELS)
    assert sorted(f.name for f in inputs.l1b_files) == sorted(_LOCAL_NAMES)
    assert inputs.emissivity_path is not None
    assert inputs.background == "semianalytic"


def test_resolve_mongolia_missing_band_raises(tmp_path):
    _make_tree(tmp_path, skip_channel="B13")
    with pytest.raises(FileNotFoundError, match="fetch_case"):
        run_case.resolve_inputs(_CASE, data_dir=tmp_path)


def test_resolve_composite_on_ahi_not_implemented(tmp_path):
    with pytest.raises(NotImplementedError):
        run_case.resolve_inputs(_CASE, data_dir=tmp_path, background="composite")


def test_abi_resolution_unchanged(tmp_path):
    # The rename must not disturb the ABI path or its error semantics.
    with pytest.raises(FileNotFoundError, match="fetch_case"):
        run_case.resolve_inputs("2017-03-23-swus", data_dir=tmp_path)


# --- load_scene bbox signature ------------------------------------------------


def test_load_scene_has_bbox_keyword():
    from shachen.io.satellite import load_scene

    params = inspect.signature(load_scene).parameters
    assert "bbox" in params
    assert params["bbox"].default is None
    assert "roles" in params  # composite-background interface intact


# --- data-gated tests on the real fetched case --------------------------------

_DATA = _ROOT / "data"
_AHI_DIR = _DATA / _CASE / "ahi"
_MERRA16 = _DATA / "merra2" / "merra2_ts_20160421.nc"
_EMIS16 = _DATA / "emissivity" / "CAM5K30EM_201604.nc"

_HAVE_SCENE = _AHI_DIR.exists() and len(list(_AHI_DIR.glob("*.DAT*"))) >= 8
_HAVE_CASE = _HAVE_SCENE and _MERRA16.exists() and _EMIS16.exists()


@pytest.mark.skipif(
    not _HAVE_SCENE,
    reason="2016-04-21 Mongolia AHI files not downloaded "
    "(run scripts/fetch_case.py 2016-04-21-mongolia)",
)
def test_ahi_scene_loads_cropped():
    from shachen.io.satellite import load_scene

    ds = load_scene(sorted(_AHI_DIR.glob("*.DAT*")), reader="ahi_hsd", bbox=_BBOX)
    assert sorted(ds.data_vars) == [
        "bt_swir_39",
        "bt_tir_104",
        "bt_tir_123",
        "bt_tir_86",
        "bt_wv_62",
        "refl_nir_160",
        "refl_vis_064",
    ]
    shapes = {v.shape for v in ds.data_vars.values()}
    assert len(shapes) == 1
    area = ds.attrs["area"]
    assert area.pixel_size_x == pytest.approx(2000, abs=10)

    # Scene.crop(ll_bbox=...) keeps the projection-space bounding rectangle
    # of the lat/lon box, so pixels outside the box survive at the corners;
    # the contract is coverage of the box plus a hard shrink versus the
    # 5500 x 5500 full disk.
    lons, lats = area.get_lonlats()
    ok = np.isfinite(lons) & np.isfinite(lats)
    assert ok.any()
    assert lons[ok].min() <= _BBOX[0] and lons[ok].max() >= _BBOX[2]
    assert lats[ok].min() <= _BBOX[1] and lats[ok].max() >= _BBOX[3]
    ny, nx = ds["bt_tir_104"].shape
    assert ny * nx < 5500 * 5500 / 4

    bt = ds["bt_tir_104"]
    assert float(bt.min()) > 180.0 and float(bt.max()) < 340.0
    assert ds.attrs["start_time"].replace(second=0, microsecond=0) == _WHEN


@pytest.mark.skipif(
    not _HAVE_CASE,
    reason="2016-04-21 Mongolia case incomplete (need AHI + MERRA-2 + CAMEL; "
    "run scripts/fetch_case.py 2016-04-21-mongolia)",
)
def test_mongolia_end_to_end(tmp_path):
    import xarray as xr

    run_case.main([_CASE, "--outdir", str(tmp_path), "--tuning", "paper"])
    ncs = list(tmp_path.glob("*.nc"))
    pngs = list(tmp_path.glob("*.png"))
    # The DEBRA composite plus the classic Dust RGB baseline PNG.
    assert len(ncs) == 1 and len(pngs) == 2
    for png in pngs:
        assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"

    out = xr.open_dataset(ncs[0])
    assert "cf_comb" in out and "rgb" in out
    cf = np.asarray(out["cf_comb"].values, dtype=float)
    assert np.isfinite(cf).any()
    rgb = np.asarray(out["rgb"].values, dtype=float)
    finite = np.isfinite(rgb)
    assert finite.any()
    assert ((rgb[finite] >= 0.0) & (rgb[finite] <= 1.0)).all()

    # Quantitative acceptance, fixed after visual validation against paper
    # Fig. 5: the plume core over the Mongolia/Inner
    # Mongolia border lights up (measured mean 0.147, p90 0.485, 40% of
    # pixels > 0.1, red-blue 0.11) while clear ocean stays dark (mean 0.005).
    from shachen.io.satellite import load_scene

    scene = load_scene(sorted(_AHI_DIR.glob("*.DAT*")), reader="ahi_hsd", bbox=_BBOX)
    lons, lats = scene.attrs["area"].get_lonlats()
    ok = np.isfinite(cf) & np.isfinite(lons) & np.isfinite(lats)
    core = ok & (lats >= 45.0) & (lats <= 47.0) & (lons >= 114.0) & (lons <= 122.0)
    ocean = ok & (lats >= 32.0) & (lats <= 36.0) & (lons >= 128.0) & (lons <= 133.0)
    assert core.any() and ocean.any()
    assert cf[core].mean() > 0.10
    assert np.percentile(cf[core], 90) > 0.35
    assert (cf[core] > 0.1).mean() > 0.25
    assert cf[ocean].mean() < 0.02

    red, blu = rgb[..., 0], rgb[..., 2]
    core_rgb = core & np.isfinite(red) & np.isfinite(blu)
    assert (red[core_rgb] - blu[core_rgb]).mean() > 0.05  # yellow dust


def test_ahi_segments_for_bbox_picks_the_northern_strips():
    """The north-China dust domain lives in the top three strips of the disk."""
    assert fetch_case.ahi_segments_for_bbox((80.0, 30.0, 135.0, 55.0)) == (1, 2, 3)


def test_ahi_segments_for_bbox_follows_the_disk_south():
    """Strip numbers run north to south, and None means the whole disk."""
    assert fetch_case.ahi_segments_for_bbox((130.0, -5.0, 150.0, 5.0)) == (5, 6)
    assert fetch_case.ahi_segments_for_bbox((110.0, -40.0, 155.0, -10.0)) == (6, 7, 8, 9)
    assert fetch_case.ahi_segments_for_bbox(None) == tuple(range(1, fetch_case.AHI_SEGMENTS + 1))


@pytest.mark.parametrize(
    "bbox",
    [
        (80.0, 30.0, 135.0, 55.0),  # the north-China dust domain
        (138.0, 30.0, 143.0, 55.0),  # the same latitudes under the satellite
        (130.0, -5.0, 150.0, 5.0),  # across the equator
        (110.0, -40.0, 155.0, -10.0),  # southern hemisphere
        (120.0, 20.0, 125.0, 25.0),  # small, mid-disk
    ],
)
def test_ahi_segments_for_bbox_covers_the_domain(bbox):
    """Whatever strips come back must span the domain's projected y range.

    Computed here from the projection rather than from latitude, which is the
    property that matters: a strip short at either end silently truncates the
    crop, and no test of specific strip numbers would catch it everywhere.
    """
    from pyproj import CRS, Transformer

    transformer = Transformer.from_crs(
        CRS.from_epsg(4326), CRS.from_dict(fetch_case.AHI_PROJ), always_xy=True
    )
    lon_min, lat_min, lon_max, lat_max = bbox
    ys = [
        transformer.transform(lon, lat)[1]
        for lon in np.linspace(lon_min, lon_max, 25)
        for lat in (lat_min, lat_max)
    ]

    segments = fetch_case.ahi_segments_for_bbox(bbox)
    strip = 2 * fetch_case.AHI_FLDK_Y_EXTENT / fetch_case.AHI_SEGMENTS
    top = fetch_case.AHI_FLDK_Y_EXTENT - (segments[0] - 1) * strip
    bottom = fetch_case.AHI_FLDK_Y_EXTENT - segments[-1] * strip
    assert bottom <= min(ys) and max(ys) <= top
    assert list(segments) == list(range(segments[0], segments[-1] + 1))  # contiguous


def test_select_ahi_keys_keeps_only_the_requested_segments():
    names = tuple(f"HS_H08_20210315_0400_B13_FLDK_R20_S{k:02d}10.DAT.bz2" for k in range(1, 11))
    selected = fetch_case.select_ahi_keys(_listing(names), channels=("B13",), segments=(1, 2, 3))
    assert [Path(k).name[-14:-8] for k in selected] == ["_S0110", "_S0210", "_S0310"]


def test_select_ahi_keys_never_drops_a_whole_disk_file():
    """``_S0101`` is one file for the entire disk, not the northern strip.

    The early-years AWS repackaging ships a timeline that way; filtering it
    down to the northern strips would silently discard the whole scene.
    """
    names = ("HS_H08_20160421_0800_B13_FLDK_R20_S0101.DAT.bz2",)
    selected = fetch_case.select_ahi_keys(_listing(names), channels=("B13",), segments=(1, 2, 3))
    assert selected == _listing(names)
    southern = fetch_case.select_ahi_keys(_listing(names), channels=("B13",), segments=(7, 8))
    assert southern == _listing(names)


def test_select_ahi_keys_raises_when_the_segments_hold_nothing():
    names = tuple(f"HS_H08_20210315_0400_B13_FLDK_R20_S{k:02d}10.DAT.bz2" for k in range(7, 11))
    with pytest.raises(RuntimeError, match="B13"):
        fetch_case.select_ahi_keys(_listing(names), channels=("B13",), segments=(1, 2))
