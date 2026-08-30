"""Fetch all inputs for a named DEBRA case: L1b bands (ABI netCDF or AHI HSD,
anonymous S3), MERRA-2 skin temperature, and CAMEL emissivity (both via
Earthdata netrc).

Usage:
    uv run python scripts/fetch_case.py 2017-03-23-swus [--quicklook] [--no-ancillary]

The --quicklook flag renders the BT 10.3/10.4 um scene to
data/<case>/quicklook_bt104.png.
"""

import argparse
import datetime as dt
import re
from dataclasses import dataclass
from pathlib import Path

from shachen.constants import COMPOSITE_WINDOW_DAYS

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

#: The 7 DEBRA roles + C14/B14 (11.2 um, classic Dust RGB green gun only),
#: constants.ABI_BANDS/AHI_BANDS order.
ABI_CHANNELS = ("C02", "C05", "C07", "C08", "C11", "C13", "C14", "C15")

AHI_CHANNELS = ("B03", "B05", "B07", "B08", "B11", "B13", "B14", "B15")

#: ABI channel names of shachen.composite.COMPOSITE_BANDS (8.4/10.3/12.3 um) -
#: the only bands a composite background day needs.
COMPOSITE_CHANNELS = ("C11", "C13", "C15")

#: Sensor dispatch tables: satpy reader, on-disk L1b glob, and the
#: per-sensor channel tuple. run_case.py imports these too.
SENSOR_READERS: dict[str, str] = {"abi": "abi_l1b", "ahi": "ahi_hsd"}
SENSOR_GLOBS: dict[str, str] = {"abi": "*.nc", "ahi": "*.DAT*"}
SENSOR_CHANNELS: dict[str, tuple[str, ...]] = {
    "abi": ABI_CHANNELS,
    "ahi": AHI_CHANNELS,
}


@dataclass(frozen=True)
class Case:
    """A reference dust case (see plan §2.4).

    ``sensor`` selects the fetch path and satpy reader (``abi`` or ``ahi``);
    ``bbox`` = (lon_min, lat_min, lon_max, lat_max) crops the loaded scene
    (None = full scene) — used to trim AHI full disks to the case domain.
    """

    bucket: str
    product: str
    when: dt.datetime  # target scan time, UTC
    sensor: str = "abi"
    bbox: tuple[float, float, float, float] | None = None


CASES: dict[str, Case] = {
    # Paper Figure 6: TX/NM dust storm, GOES-16 (Mode 3, 5-min CONUS).
    "2017-03-23-swus": Case("noaa-goes16", "ABI-L1b-RadC", dt.datetime(2017, 3, 23, 21, 45)),
    # CIRA Quick Guide reference images, Southern Plains.
    "2020-12-13-splains": Case("noaa-goes16", "ABI-L1b-RadC", dt.datetime(2020, 12, 13, 21, 20)),
    "2020-12-23-splains": Case("noaa-goes16", "ABI-L1b-RadC", dt.datetime(2020, 12, 23, 17, 41)),
    # Paper Figures 3 & 5: Mongolia -> China dust storm, Himawari-8 full disk
    # cropped to the papers' display domain.
    "2016-04-21-mongolia": Case(
        "noaa-himawari8",
        "AHI-L1b-FLDK",
        dt.datetime(2016, 4, 21, 8, 0),
        sensor="ahi",
        bbox=(80.0, 30.0, 135.0, 55.0),
    ),
}


def _scan_start(key: str) -> dt.datetime:
    """Parse the sYYYYJJJHHMMSSs field of an ABI filename."""
    stamp = key.split("_s")[1][:13]
    return dt.datetime.strptime(stamp, "%Y%j%H%M%S")


def fetch_abi(case: Case, out_dir: Path) -> list[Path]:
    """Download the ``ABI_CHANNELS`` bands for the scan nearest case.when."""
    import s3fs

    fs = s3fs.S3FileSystem(anon=True)
    prefix = f"{case.bucket}/{case.product}/{case.when:%Y/%j/%H}/"
    keys = fs.ls(prefix)

    out_dir.mkdir(parents=True, exist_ok=True)
    local: list[Path] = []
    for chan in ABI_CHANNELS:
        candidates = [k for k in keys if f"{chan}_" in k]
        if not candidates:
            raise RuntimeError(f"No {chan} files under s3://{prefix}")
        best = min(candidates, key=lambda k: abs(_scan_start(k) - case.when))
        dest = out_dir / Path(best).name
        if not dest.exists():
            print(f"  {chan}: s3://{best}")
            fs.get(best, str(dest))
        local.append(dest)
    return local


def ahi_timeline_prefix(case: Case) -> str:
    """S3 prefix of the 10-minute AHI full-disk timeline holding ``case.when``.

    ``f"{bucket}/{product}/%Y/%m/%d/%H%M/"`` — the AWS Himawari layout has
    one directory per timeline. ``case.when`` must sit exactly on a timeline
    (minute a multiple of 10, second zero); anything else raises ValueError.
    """
    when = case.when
    if when.minute % 10 != 0 or when.second != 0 or when.microsecond != 0:
        raise ValueError(f"{when:%Y-%m-%d %H:%M:%S} is not on a 10-minute AHI timeline")
    return f"{case.bucket}/{case.product}/{when:%Y/%m/%d/%H%M}/"


_AHI_RESOLUTION_RE = re.compile(r"_R(\d{2})_")

#: ``_S{segment:02d}{total:02d}`` in an HSD basename: which strip of how many.
#: The field ends the stem, so it is followed by the extension, not by "_".
_AHI_SEGMENT_RE = re.compile(r"_S(\d{2})(\d{2})(?=[._])")

#: A full-disk HSD timeline is cut into 10 strips, S01 northernmost, each a
#: tenth of the 2-km fixed grid's 5500 lines (+-5,500,000 m in y).
AHI_SEGMENTS = 10
AHI_FLDK_Y_EXTENT = 5_500_000.0

#: The AHI fixed grid's projection (Himawari-8/9 parked at 140.7E).
AHI_PROJ = {
    "proj": "geos",
    "h": 35785863.0,
    "lon_0": 140.7,
    "sweep": "y",
    "a": 6378137.0,
    "rf": 298.257024882273,
}

#: Slop added to a domain's projected extent before it is turned into strip
#: numbers, so an ellipsoid or pixel-center difference can never truncate the
#: crop. One strip is 1,100,000 m, so this is cheap insurance.
_SEGMENT_MARGIN_M = 20_000.0


def ahi_segments_for_bbox(bbox: tuple[float, float, float, float] | None) -> tuple[int, ...]:
    """The 1-based HSD strip numbers whose data covers ``bbox`` (None -> all).

    A full disk arrives as 10 north-to-south strips, and a regional domain
    needs only a few: the north-China dust domain lives in S01-S03, so
    fetching all ten downloads about four times what the crop keeps.

    The domain's boundary is projected onto the fixed grid rather than read
    off latitude alone — in a geostationary projection y depends on longitude
    too, so a domain's northernmost *projected* point need not sit at its
    northernmost corner.
    """
    if bbox is None:
        return tuple(range(1, AHI_SEGMENTS + 1))
    from pyproj import CRS, Transformer

    lon_min, lat_min, lon_max, lat_max = bbox
    steps = 100
    lons = [lon_min + (lon_max - lon_min) * i / steps for i in range(steps + 1)]
    lats = [lat_min + (lat_max - lat_min) * i / steps for i in range(steps + 1)]
    boundary = (
        [(lon, lat_min) for lon in lons]
        + [(lon, lat_max) for lon in lons]
        + [(lon_min, lat) for lat in lats]
        + [(lon_max, lat) for lat in lats]
    )
    transformer = Transformer.from_crs(CRS.from_epsg(4326), CRS.from_dict(AHI_PROJ), always_xy=True)
    ys = [
        y
        for _, y in (transformer.transform(lon, lat) for lon, lat in boundary)
        if abs(y) != float("inf") and y == y
    ]
    if not ys:
        raise RuntimeError(f"bbox {bbox} does not intersect the AHI disk")

    strip = 2 * AHI_FLDK_Y_EXTENT / AHI_SEGMENTS
    first = int((AHI_FLDK_Y_EXTENT - (max(ys) + _SEGMENT_MARGIN_M)) // strip) + 1
    last = int((AHI_FLDK_Y_EXTENT - (min(ys) - _SEGMENT_MARGIN_M)) // strip) + 1
    return tuple(range(max(1, first), min(AHI_SEGMENTS, last) + 1))


def _in_segments(key: str, segments: tuple[int, ...]) -> bool:
    """Whether an HSD key belongs to ``segments`` — single-file disks always do.

    ``_S0110_`` is strip 1 of 10; ``_S0101_`` is the whole disk in one file
    (the early-years AWS repackaging). Only a genuinely segmented timeline is
    filtered, so restricting strips can never drop a whole-disk file.
    """
    match = _AHI_SEGMENT_RE.search(Path(key).name)
    if match is None:
        return True
    segment, total = int(match.group(1)), int(match.group(2))
    return total <= 1 or segment in segments


def select_ahi_keys(
    keys: list[str],
    channels: tuple[str, ...] = AHI_CHANNELS,
    segments: tuple[int, ...] | None = None,
) -> list[str]:
    """Pick each channel's HSD keys for ``channels`` from a timeline listing.

    For each channel, in order: every key containing ``_{chan}_`` and
    ``_R20_`` (sorted; the early-years AWS repackaging is single-segment but
    multiple segments are tolerated). Timelines from ~2019 on carry native
    HSD only — no 2-km repackaging, B03 at R05 / B05 at R10 in 10 segments —
    so with no R20 match the coarsest resolution present for that channel is
    taken instead (load_scene's native resample aggregates it onto the 2-km
    grid). A channel with no keys at all raises RuntimeError naming it.

    ``segments`` (from :func:`ahi_segments_for_bbox`) keeps only the strips a
    regional domain needs; None keeps the whole disk. A channel whose keys are
    all outside the requested strips still raises RuntimeError — a domain the
    timeline does not cover is a mistake worth hearing about.
    """
    selected: list[str] = []
    for chan in channels:
        candidates = [k for k in keys if f"_{chan}_" in k]
        if segments is not None:
            candidates = [k for k in candidates if _in_segments(k, segments)]
        matches = sorted(k for k in candidates if "_R20_" in k)
        if not matches and candidates:
            coarsest = max(
                match.group(1)
                for k in candidates
                if (match := _AHI_RESOLUTION_RE.search(Path(k).name))
            )
            matches = sorted(k for k in candidates if f"_R{coarsest}_" in k)
        if not matches:
            raise RuntimeError(f"No {chan} files in timeline listing")
        selected.extend(matches)
    return selected


def ahi_local_name(key: str) -> str:
    """Local filename an AHI S3 key is stored under.

    satpy's ahi_hsd ``.DAT.bz2`` file patterns hardcode each band's *native*
    resolution, so the repackaged 2-km B03 (natively 0.5 km) is unreadable
    compressed — only the uncompressed ``.DAT`` pattern is
    resolution-generic. Repackaged R20 B03 keys therefore lose their
    ``.bz2`` suffix (:func:`fetch_ahi` stores them decompressed); native-
    resolution B03 (R05) matches the compressed pattern as-is, so it — like
    every other band — keeps the downloaded basename.
    """
    name = Path(key).name
    if "_B03_" in name and "_R20_" in name and name.endswith(".bz2"):
        return name[: -len(".bz2")]
    return name


def fetch_ahi(
    case: Case,
    out_dir: Path,
    channels: tuple[str, ...] = AHI_CHANNELS,
    segments: tuple[int, ...] | None = None,
) -> list[Path]:
    """Download ``channels`` (2-km R20 files) of ``case``'s timeline.

    Lists :func:`ahi_timeline_prefix` (missing/empty prefix -> RuntimeError
    naming it), selects ``channels`` via :func:`select_ahi_keys` (a caller that
    only feeds DEBRA can pin the 7 bands here), downloads files not
    already in ``out_dir``, and returns the local paths in selection order.
    Each key lands under :func:`ahi_local_name` — B03 is decompressed on
    arrival (stdlib bz2, the compressed download removed) and the
    exists-check runs against that final name, so re-runs stay idempotent.

    ``segments`` defaults to the strips ``case.bbox`` needs
    (:func:`ahi_segments_for_bbox`), which is what the crop keeps anyway —
    a quarter of the bytes for a regional domain. Pass a tuple to override,
    or ``range(1, AHI_SEGMENTS + 1)`` for the whole disk. **The cache
    directory therefore holds one bbox's strips**: widening a domain later
    means refetching, since the exists-check cannot know what is missing.
    """
    import bz2

    import s3fs

    fs = s3fs.S3FileSystem(anon=True)
    prefix = ahi_timeline_prefix(case)
    try:
        keys = fs.ls(prefix)
    except FileNotFoundError:
        keys = []
    if not keys:
        raise RuntimeError(f"Nothing under s3://{prefix}")

    if segments is None:
        segments = ahi_segments_for_bbox(case.bbox)
    out_dir.mkdir(parents=True, exist_ok=True)
    local: list[Path] = []
    for key in select_ahi_keys(keys, channels, segments):
        chan = Path(key).name.split("_")[4]
        dest = out_dir / ahi_local_name(key)
        if not dest.exists():
            print(f"  {chan}: s3://{key}")
            # Download to a .part name and rename once complete: an
            # interrupted transfer must never leave a truncated file under
            # the final name, or the exists-check would trust it on resume.
            download = out_dir / (Path(key).name + ".part")
            fs.get(key, str(download))
            if ahi_local_name(key) != Path(key).name:
                # B03: satpy's bz2 patterns can't read the 2-km repackaging,
                # so store it decompressed and drop the compressed download.
                decompressed = dest.with_suffix(dest.suffix + ".part")
                decompressed.write_bytes(bz2.decompress(download.read_bytes()))
                decompressed.replace(dest)
                download.unlink()
            else:
                download.replace(dest)
        local.append(dest)
    return local


def fetch_l1b(case: Case, out_dir: Path) -> list[Path]:
    """Fetch ``case``'s L1b files, dispatching on ``case.sensor``.

    ``abi`` -> :func:`fetch_abi`, ``ahi`` -> :func:`fetch_ahi`; any other
    sensor raises ValueError.
    """
    # Looked up at call time (module globals) so tests can monkeypatch them.
    fetchers = {"abi": fetch_abi, "ahi": fetch_ahi}
    if case.sensor not in fetchers:
        raise ValueError(f"unknown sensor {case.sensor!r} (choose from {sorted(fetchers)})")
    return fetchers[case.sensor](case, out_dir)


def composite_days(when: dt.datetime, days: int) -> list[dt.date]:
    """The ``days`` calendar days preceding ``when.date()``, chronological.

    The case day itself is excluded (a dusty day can be the warmest pixel and
    poison the composite): ``composite_days(2020-12-13 21:20, 14)`` is
    2020-11-29 .. 2020-12-12.
    """
    return [when.date() - dt.timedelta(days=n) for n in range(days, 0, -1)]


def fetch_composite(case: Case, out_dir: Path, days: int = COMPOSITE_WINDOW_DAYS) -> list[Path]:
    """Download the TIR composite window for ``case``.

    For each day of :func:`composite_days`, fetch the ``COMPOSITE_CHANNELS``
    files of the scan nearest the case's time of day into
    ``out_dir/<YYYYMMDD>/`` (files already on disk are kept). A day with
    nothing on S3 is skipped with a printed warning. Returns the per-day
    directories actually populated, chronological.
    """
    import s3fs

    fs = s3fs.S3FileSystem(anon=True)
    populated: list[Path] = []
    for day in composite_days(case.when, days):
        target = dt.datetime.combine(day, case.when.time())
        prefix = f"{case.bucket}/{case.product}/{target:%Y/%j/%H}/"
        try:
            keys = fs.ls(prefix)
        except FileNotFoundError:
            keys = []
        if not keys:
            print(f"  {day:%Y%m%d}: nothing under s3://{prefix}, skipping")
            continue

        day_dir = out_dir / f"{day:%Y%m%d}"
        day_dir.mkdir(parents=True, exist_ok=True)
        for chan in COMPOSITE_CHANNELS:
            candidates = [k for k in keys if f"{chan}_" in k]
            if not candidates:
                raise RuntimeError(f"No {chan} files under s3://{prefix}")
            best = min(candidates, key=lambda k: abs(_scan_start(k) - target))
            dest = day_dir / Path(best).name
            if not dest.exists():
                print(f"  {day:%Y%m%d} {chan}: s3://{best}")
                fs.get(best, str(dest))
        populated.append(day_dir)
    return populated


def quicklook(
    files: list[Path],
    out_png: Path,
    reader: str = "abi_l1b",
    bbox: tuple[float, float, float, float] | None = None,
) -> None:
    """Render the BT 10.3/10.4 um field with coastlines for geographic sanity."""
    import cartopy.feature as cfeature
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from shachen.constants import Band
    from shachen.io.satellite import load_scene

    ds = load_scene(files, reader=reader, bbox=bbox)
    bt = ds[f"bt_{Band.TIR_104.value}"]
    crs = ds.attrs["area"].to_cartopy_crs()

    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={"projection": crs}, layout="constrained")
    im = ax.imshow(
        bt.values,
        extent=crs.bounds,
        origin="upper",
        cmap="gray_r",
        vmin=220,
        vmax=310,
        transform=crs,
    )
    ax.coastlines(resolution="50m", color="tab:red", linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, edgecolor="tab:red", linewidth=0.6)
    ax.add_feature(cfeature.STATES, edgecolor="tab:red", linewidth=0.3)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="tab:blue")
    gl.top_labels = gl.right_labels = False
    fig.colorbar(im, ax=ax, shrink=0.7, label="BT 10.3/10.4 µm [K]")
    ax.set_title(f"BT 10.3/10.4 µm — {ds.attrs['start_time']:%Y-%m-%d %H:%M} UTC")
    fig.savefig(out_png, dpi=150)
    print(f"  quicklook -> {out_png}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=sorted(CASES), nargs="?", default="2017-03-23-swus")
    parser.add_argument("--quicklook", action="store_true", help="render BT10.3 PNG")
    parser.add_argument("--no-ancillary", action="store_true", help="skip MERRA-2 + emissivity")
    parser.add_argument(
        "--composite",
        type=int,
        default=0,
        metavar="DAYS",
        help="also fetch DAYS preceding days of TIR bands for the "
        "composite background (0 = off; default is "
        f"{COMPOSITE_WINDOW_DAYS})",
    )
    args = parser.parse_args()

    case = CASES[args.case]
    case_dir = DATA_DIR / args.case

    if args.composite and case.sensor != "abi":
        parser.error("--composite is ABI-only")

    print(f"[1/3] {case.sensor.upper()} L1b bands for {case.when:%Y-%m-%d %H:%M} UTC:")
    files = fetch_l1b(case, case_dir / case.sensor)

    if args.composite:
        print(f"[+] composite window ({args.composite} days, TIR only):")
        fetched = fetch_composite(case, case_dir / "composite", args.composite)
        print(f"  {len(fetched)} day(s) -> {case_dir / 'composite'}")

    if not args.no_ancillary:
        from shachen.io.emissivity import fetch_emissivity
        from shachen.io.merra import fetch_skin_temperature

        print("[2/3] MERRA-2 skin temperature (Earthdata netrc login):")
        ts_path = fetch_skin_temperature(case.when.date(), DATA_DIR / "merra2")
        print(f"  TS cache -> {ts_path}")

        print("[3/3] CAMEL monthly emissivity:")
        emis_path = fetch_emissivity(case.when.date(), DATA_DIR / "emissivity")
        print(f"  emissivity -> {emis_path}")

    if args.quicklook:
        quicklook(
            files,
            case_dir / "quicklook_bt104.png",
            reader=SENSOR_READERS[case.sensor],
            bbox=case.bbox,
        )


if __name__ == "__main__":
    main()
