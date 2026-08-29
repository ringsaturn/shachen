"""Cloud-cleared composite clear-sky background (Miller et al. 2017, §3.2).

The emissivity-free alternative to :mod:`shachen.background`: stack the same
time-of-day scenes from the preceding ~14 days and, per pixel, keep the three
TIR brightness temperatures from the day with the warmest BT10.4 (clouds are
cold; the warmest day is the best clear-sky estimate). Because the composite
is built from real observations it carries the split-window water-vapor
depression the semianalytic (emissivity x Planck, no atmosphere) background
lacks — the ~1 K high bias that zeroes DT1/DT2 on transparent winter plumes
(see docs/deviations.md).

Spectral coherence: all three bands come from the *same* selected day, and a
day is a candidate at a pixel only where all three bands are finite there.
"""

from collections.abc import Sequence

import numpy as np
import xarray as xr

from shachen.constants import Band

#: The bands a composite scene must carry (the three DEBRA TIR windows).
COMPOSITE_BANDS: tuple[Band, ...] = (Band.TIR_86, Band.TIR_104, Band.TIR_123)


def composite_background(scenes: Sequence[xr.Dataset]) -> xr.Dataset:
    """Clear-sky background signals from a stack of same-time-of-day scenes.

    Each element of ``scenes`` must carry 2-D ``bt_tir_86``, ``bt_tir_104``,
    ``bt_tir_123`` (K) on one common grid; ValueError on an empty sequence, a
    missing variable, or any shape mismatch. Per pixel, a scene is a
    *candidate* only where all three bands are finite; the candidate with the
    warmest ``bt_tir_104`` is selected (first wins on exact ties) and all
    three bands are taken from it. Pixels with zero candidates are NaN.

    Returns the :func:`shachen.background.background_signals` contract —
    ``bt_bg_tir_86``, ``bt_bg_tir_104``, ``bt_bg_tir_123`` (units "K"),
    ``rsw_bg`` (= bt_bg_tir_123 - bt_bg_tir_104) and ``btd_bg``
    (= bt_bg_tir_86 - bt_bg_tir_104) — plus ``n_valid`` (per-pixel candidate
    count, integer dtype) and ``attrs["n_scenes"] = len(scenes)``.
    """
    if len(scenes) == 0:
        raise ValueError("composite_background needs at least one scene")

    # Gather each band as a (n_scenes, y, x) stack, validating each.
    shape: tuple[int, ...] | None = None
    stacks: dict[Band, np.ndarray] = {}
    for band in COMPOSITE_BANDS:
        name = f"bt_{band.value}"
        layers = []
        for i, scene in enumerate(scenes):
            if name not in scene:
                raise ValueError(f"scene {i} is missing variable {name!r}")
            values = np.asarray(scene[name].values, dtype=float)
            if shape is None:
                shape = values.shape
            elif values.shape != shape:
                raise ValueError(
                    f"shape mismatch: scene {i} {name} {values.shape} vs expected {shape}"
                )
            layers.append(values)
        stacks[band] = np.stack(layers, axis=0)

    # A scene is a candidate at a pixel only where all three bands are finite.
    candidate = np.logical_and.reduce([np.isfinite(stacks[band]) for band in COMPOSITE_BANDS])
    n_valid = candidate.sum(axis=0)
    any_candidate = n_valid > 0

    # Warmest-BT10.4 candidate; argmax takes the first scene on exact ties.
    # Non-candidates are masked to -inf so they can never win; all-masked
    # columns fall back to index 0 and are NaN-ed via any_candidate below.
    bt104_masked = np.where(candidate, stacks[Band.TIR_104], -np.inf)
    selected = bt104_masked.argmax(axis=0)

    out = xr.Dataset(attrs={"n_scenes": len(scenes)})
    for band in COMPOSITE_BANDS:
        picked = np.take_along_axis(stacks[band], selected[np.newaxis], axis=0)[0]
        out[f"bt_bg_{band.value}"] = (
            ("y", "x"),
            np.where(any_candidate, picked, np.nan),
        )
    out["rsw_bg"] = out["bt_bg_tir_123"] - out["bt_bg_tir_104"]
    out["btd_bg"] = out["bt_bg_tir_86"] - out["bt_bg_tir_104"]
    out["n_valid"] = (("y", "x"), n_valid)

    for band in COMPOSITE_BANDS:
        out[f"bt_bg_{band.value}"].attrs.update(
            units="K",
            long_name=f"clear-sky background brightness temperature ({band.value})",
        )
    out["rsw_bg"].attrs.update(
        units="K", long_name="background split-window signal BT(12.3) - BT(10.4)"
    )
    out["btd_bg"].attrs.update(
        units="K", long_name="background brightness temperature difference BT(8.6) - BT(10.4)"
    )
    out["n_valid"].attrs.update(long_name="number of candidate (all-bands-finite) scenes per pixel")
    return out
