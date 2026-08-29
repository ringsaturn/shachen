"""Render the DEBRA RGB composite to a georeferenced PNG with map overlays.

Paper section 4.2: "Georeferenced metainformation such as coastlines,
political boundaries, latitude/longitude grid, and date/time information ...
are overlaid to the R/G/B imagery." Uses cartopy on the scene's projection,
same styling as the quicklook renderer.
"""

import datetime as dt
from pathlib import Path

import numpy as np
from pyresample.geometry import AreaDefinition


def render_debra_png(
    rgb: np.ndarray,
    area: AreaDefinition,
    when: dt.datetime,
    out_path: Path,
    title: str | None = None,
    dpi: int = 150,
) -> Path:
    """Write ``rgb`` as a PNG with coastlines, borders, and a lat/lon grid.

    ``rgb`` must be a ``(y, x, 3)`` uint8 array (from
    :func:`shachen.imagery.to_uint8`) on ``area``'s grid — anything else raises
    ValueError. ``when`` is the scan start time (UTC), used in the default
    title ``"DEBRA-Dust — %Y-%m-%d %H:%M UTC"``. Renders on the Agg backend
    and returns ``out_path``.
    """
    _validate_rgb(rgb)

    import cartopy.feature as cfeature
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    crs = area.to_cartopy_crs()

    fig, ax = plt.subplots(figsize=_FIGSIZE, subplot_kw={"projection": crs}, layout="constrained")
    try:
        ax.imshow(rgb, extent=crs.bounds, origin="upper", transform=crs)
        # Paper section 4.2: coastlines, political boundaries, lat/lon grid,
        # date/time information overlaid on the R/G/B imagery.
        ax.coastlines(resolution="50m", color="tab:red", linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, edgecolor="tab:red", linewidth=0.6)
        ax.add_feature(cfeature.STATES, edgecolor="tab:red", linewidth=0.3)
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="tab:blue")
        gl.top_labels = gl.right_labels = False
        ax.set_title(title if title is not None else _default_title(when))
        fig.savefig(out_path, dpi=dpi)
    finally:
        plt.close(fig)
    return out_path


_FIGSIZE = (12.0, 8.0)


def _default_title(when: dt.datetime) -> str:
    """Date/time annotation required by paper section 4.2."""
    return f"DEBRA-Dust — {when:%Y-%m-%d %H:%M} UTC"


def _validate_rgb(rgb: np.ndarray) -> None:
    """Require a ``(y, x, 3)`` uint8 array, as produced by ``imagery.to_uint8``."""
    arr = np.asarray(rgb)
    if arr.ndim != 3 or arr.shape[-1] != 3:
        raise ValueError(f"rgb must have shape (y, x, 3), got {arr.shape}")
    if arr.dtype != np.uint8:
        raise ValueError(f"rgb must be uint8, got dtype {arr.dtype}")
