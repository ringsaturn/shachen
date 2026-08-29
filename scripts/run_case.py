"""Run DEBRA end-to-end on a named reference case -> netCDF + PNG.

Usage:
    uv run python scripts/run_case.py 2017-03-23-swus [--color yellow]
        [--tuning abi] [--outdir data/2017-03-23-swus] [--dpi 150]

Inputs must already be on disk (run scripts/fetch_case.py first); outputs are
a netCDF carrying CF_comb, all intermediate fields, and the float RGB
composite, plus the rendered PNG.
"""

import argparse
import datetime as dt
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from fetch_case import (  # noqa: E402
    CASES,
    COMPOSITE_CHANNELS,
    DATA_DIR,
    SENSOR_CHANNELS,
    SENSOR_GLOBS,
    SENSOR_READERS,
)

from shachen.constants import (  # noqa: E402
    ABI_TUNED,
    COLOR_DIMMING,
    COMPOSITE_MIN_DAYS,
    COMPOSITE_WINDOW_DAYS,
    DEFAULTS,
    DebraConstants,
)

#: Selectable constant sets. ``abi`` is the ABI retune (the default);
#: ``paper`` reverts to the published Miller et al. 2017/2020 values.
TUNING: dict[str, DebraConstants] = {"abi": ABI_TUNED, "paper": DEFAULTS}

#: Selectable background modes: ``semianalytic`` is the CAMEL-emissivity x
#: Planck background, ``composite`` the cloud-cleared multi-day composite
#: (needs fetch_case.py --composite first).
BACKGROUNDS: tuple[str, ...] = ("semianalytic", "composite")


@dataclass(frozen=True)
class CaseInputs:
    """Resolved on-disk inputs for one reference case.

    ``emissivity_path`` is None in composite mode (CAMEL is not needed);
    ``composite_dirs`` are the per-day TIR directories, chronological, and
    empty in semianalytic mode. ``sensor`` and ``bbox`` are copied from the
    Case: they pick the satpy reader and the lat/lon crop that
    ``run_case`` passes to ``load_scene``.
    """

    name: str
    when: dt.datetime
    l1b_files: list[Path]
    merra_path: Path
    emissivity_path: Path | None
    background: str = "semianalytic"
    composite_dirs: tuple[Path, ...] = ()
    sensor: str = "abi"
    bbox: tuple[float, float, float, float] | None = None


#: Attribute value types netCDF can store as-is.
_NC_SCALARS = (str, int, float, complex, bool, np.number, np.ndarray, list, tuple)


def _sanitize_attr(value: object) -> object:
    """Coerce one attr value into something the netCDF writer accepts."""
    if isinstance(value, (dt.datetime, dt.date)):
        return value.isoformat()
    if isinstance(value, _NC_SCALARS):
        return value
    # AreaDefinition and anything else exotic: keep the description only.
    return str(value)


def _sanitize_attrs(attrs: dict) -> dict:
    return {str(key): _sanitize_attr(value) for key, value in attrs.items()}


def resolve_inputs(
    case_name: str,
    data_dir: Path = DATA_DIR,
    background: str = "semianalytic",
) -> CaseInputs:
    """Locate the on-disk inputs for ``case_name`` in the given background mode.

    Uses the fetch_case.py layout:
    ``<data_dir>/<case>/<sensor>/<SENSOR_GLOBS[sensor]>`` L1b files (at least
    one matching ``{chan}_`` per channel of ``SENSOR_CHANNELS[sensor]``) and
    ``<data_dir>/merra2/merra2_ts_<%Y%m%d>.nc`` in both modes;
    ``semianalytic`` additionally needs
    ``<data_dir>/emissivity/CAM5K30EM_<%Y%m>.nc``, while ``composite``
    instead needs at least ``COMPOSITE_MIN_DAYS`` day directories
    ``<data_dir>/<case>/composite/<YYYYMMDD>/`` each holding >= 3 ``*.nc``
    (and sets ``emissivity_path=None``). An unknown case name raises KeyError;
    an unknown mode raises ValueError; ``composite`` on a non-ABI sensor
    raises NotImplementedError (AHI keeps the composite background
    ABI-only); missing files raise FileNotFoundError whose message points at
    ``scripts/fetch_case.py`` (including the ``--composite`` flag in
    composite mode). The returned ``CaseInputs`` carries the case's sensor
    and bbox.
    """
    if background not in BACKGROUNDS:
        raise ValueError(f"unknown background mode {background!r} (choose from {BACKGROUNDS})")
    case = CASES[case_name]  # unknown case -> KeyError
    if background == "composite" and case.sensor != "abi":
        raise NotImplementedError(
            f"composite background is ABI-only; case {case_name!r} uses sensor {case.sensor!r}"
        )
    root = Path(data_dir)

    l1b_dir = root / case_name / case.sensor
    l1b_glob = SENSOR_GLOBS[case.sensor]
    l1b_files = sorted(l1b_dir.glob(l1b_glob))
    merra_path = root / "merra2" / f"merra2_ts_{case.when:%Y%m%d}.nc"

    missing: list[str] = []
    absent = [
        chan
        for chan in SENSOR_CHANNELS[case.sensor]
        if not any(f"{chan}_" in f.name for f in l1b_files)
    ]
    if absent:
        missing.append(f"{l1b_dir}/{l1b_glob} (no file for {', '.join(absent)})")
    if not merra_path.exists():
        missing.append(str(merra_path))

    emissivity_path: Path | None = None
    composite_dirs: tuple[Path, ...] = ()
    if background == "semianalytic":
        emissivity_path = root / "emissivity" / f"CAM5K30EM_{case.when:%Y%m}.nc"
        if not emissivity_path.exists():
            missing.append(str(emissivity_path))
        fetch_hint = f"uv run python scripts/fetch_case.py {case_name}"
    else:
        composite_root = root / case_name / "composite"
        day_dirs = (
            sorted(d for d in composite_root.iterdir() if d.is_dir())
            if composite_root.exists()
            else []
        )
        composite_dirs = tuple(
            d for d in day_dirs if len(list(d.glob("*.nc"))) >= len(COMPOSITE_CHANNELS)
        )
        if len(composite_dirs) < COMPOSITE_MIN_DAYS:
            missing.append(
                f"{composite_root}/<YYYYMMDD>/*.nc (found "
                f"{len(composite_dirs)} usable day(s), need {COMPOSITE_MIN_DAYS})"
            )
        fetch_hint = (
            f"uv run python scripts/fetch_case.py {case_name} --composite {COMPOSITE_WINDOW_DAYS}"
        )
    if missing:
        listing = "\n  ".join(missing)
        raise FileNotFoundError(
            f"case {case_name!r} is missing inputs:\n  {listing}\nrun: {fetch_hint}"
        )

    return CaseInputs(
        name=case_name,
        when=case.when,
        l1b_files=l1b_files,
        merra_path=merra_path,
        emissivity_path=emissivity_path,
        background=background,
        composite_dirs=composite_dirs,
        sensor=case.sensor,
        bbox=case.bbox,
    )


def run_case(
    inputs: CaseInputs,
    color: str = "yellow",
    constants: DebraConstants = DEFAULTS,
) -> xr.Dataset:
    """Load the case inputs and run the full DEBRA chain (Eqs. 1-29).

    Chains :func:`shachen.pipeline.run_debra` and
    :func:`shachen.imagery.debra_imagery` (``color`` selects a
    ``COLOR_DIMMING`` preset; unknown names raise KeyError) and returns their
    merged variables on the scene grid, scene attrs preserved. The scene is
    loaded with ``SENSOR_READERS[inputs.sensor]`` and cropped to
    ``inputs.bbox`` when set.

    :func:`shachen.pipeline.run_dust_rgb` runs over the same scene and its
    ``dust_rgb`` variable is merged in, so every case also carries the
    classic baseline to compare the enhancement against.

    ``inputs.background`` selects the background source: ``semianalytic``
    loads CAMEL emissivity as before; ``composite`` loads each of
    ``inputs.composite_dirs`` TIR-only (``load_scene`` with
    ``roles=COMPOSITE_BANDS``), builds
    :func:`shachen.composite.composite_background`, and passes it to
    ``run_debra(background=...)``.
    """
    from shachen.imagery import debra_imagery
    from shachen.io.merra import load_skin_temperature
    from shachen.io.satellite import load_scene
    from shachen.pipeline import run_debra, run_dust_rgb

    dimming = COLOR_DIMMING[color]  # unknown color -> KeyError

    scene = load_scene(inputs.l1b_files, reader=SENSOR_READERS[inputs.sensor], bbox=inputs.bbox)
    skin_temperature = load_skin_temperature(inputs.merra_path, scene.attrs["start_time"])

    if inputs.background == "composite":
        from shachen.composite import COMPOSITE_BANDS, composite_background

        day_scenes = [
            load_scene(sorted(day_dir.glob("*.nc")), roles=COMPOSITE_BANDS)
            for day_dir in inputs.composite_dirs
        ]
        bg = composite_background(day_scenes)
        debra_out = run_debra(scene, skin_temperature, constants=constants, background=bg)
    else:
        from shachen.io.emissivity import load_band_emissivity

        emissivity = load_band_emissivity(inputs.emissivity_path)
        debra_out = run_debra(scene, skin_temperature, emissivity, constants)
    imagery = debra_imagery(scene, debra_out, constants, dimming)

    result = xr.merge([debra_out, imagery, run_dust_rgb(scene)], combine_attrs="drop")
    result.attrs.update(scene.attrs)
    return result


def save_outputs(
    result: xr.Dataset,
    out_dir: Path,
    stem: str,
    dpi: int = 150,
) -> tuple[Path, Path]:
    """Write ``result`` to ``<stem>*.nc`` and render ``<stem>*.png``.

    The netCDF carries every data variable in ``result`` with attrs sanitized
    for serialization (``area`` replaced by its string description,
    ``start_time`` by its ISO string); the PNG is rendered via
    :func:`shachen.render.render_debra_png` from ``result["rgb"]`` converted
    with :func:`shachen.imagery.to_uint8`. When ``result`` carries
    ``dust_rgb`` (the classic baseline from
    :func:`shachen.pipeline.run_dust_rgb`), a second PNG
    ``<stem>_dustrgb.png`` is rendered next to it for side-by-side
    comparison. Returns ``(nc_path, png_path)``.
    """
    from shachen.imagery import to_uint8
    from shachen.render import render_debra_png

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    nc_path = out_dir / f"{stem}.nc"
    png_path = out_dir / f"{stem}.png"

    # Serialize a copy: the caller's attrs (AreaDefinition, datetime) are still
    # needed for the map rendering below.
    to_write = result.copy(deep=False)
    to_write.attrs = _sanitize_attrs(result.attrs)
    for name in to_write.variables:
        to_write[name].attrs = _sanitize_attrs(result[name].attrs)
    to_write.to_netcdf(nc_path)

    render_debra_png(
        to_uint8(result["rgb"]),
        result.attrs["area"],
        result.attrs["start_time"],
        png_path,
        dpi=dpi,
    )
    if "dust_rgb" in result:
        render_debra_png(
            to_uint8(result["dust_rgb"]),
            result.attrs["area"],
            result.attrs["start_time"],
            out_dir / f"{stem}_dustrgb.png",
            title=f"Dust RGB — {result.attrs['start_time']:%Y-%m-%d %H:%M} UTC",
            dpi=dpi,
        )
    return nc_path, png_path


def main(argv: list[str] | None = None) -> None:
    """CLI: ``run_case.py CASE [--color NAME] [--tuning NAME] [--outdir DIR]
    [--dpi N]``.

    ``CASE`` is one of ``fetch_case.CASES``; ``--color`` one of
    ``COLOR_DIMMING`` (default ``yellow``); ``--tuning`` one of ``TUNING``
    (default ``abi``); ``--background`` one of ``BACKGROUNDS`` (default
    ``semianalytic``; ``composite`` appends ``_composite`` to the output stem
    so both runs coexist); ``--outdir`` defaults to ``data/<case>``. Prints
    the two output paths.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=sorted(CASES), help="reference case name")
    parser.add_argument(
        "--color",
        choices=sorted(COLOR_DIMMING),
        default="yellow",
        help="dust coloration preset (default: yellow)",
    )
    parser.add_argument(
        "--tuning",
        choices=sorted(TUNING),
        default="abi",
        help=("constant set: abi = ABI retune (default), paper = published Miller et al. values"),
    )
    parser.add_argument(
        "--background",
        choices=BACKGROUNDS,
        default="semianalytic",
        help=(
            "background source: semianalytic = CAMEL emissivity x Planck "
            "(default), composite = cloud-cleared multi-day composite "
            "(fetch_case.py --composite first)"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="output directory (default: data/<case>)",
    )
    parser.add_argument("--dpi", type=int, default=150, help="PNG resolution")
    args = parser.parse_args(argv)

    out_dir = args.outdir if args.outdir is not None else DATA_DIR / args.case

    inputs = resolve_inputs(args.case, background=args.background)
    result = run_case(inputs, color=args.color, constants=TUNING[args.tuning])
    stem = f"{args.case}_{inputs.when:%Y%m%dT%H%M}_debra"
    if args.background == "composite":
        stem += "_composite"
    nc_path, png_path = save_outputs(result, out_dir, stem, dpi=args.dpi)

    print(f"netCDF -> {nc_path}")
    print(f"PNG    -> {png_path}")
    print(f"PNG    -> {png_path.with_name(f'{stem}_dustrgb.png')} (classic Dust RGB baseline)")


if __name__ == "__main__":
    main()
