"""Acceptance tests for shachen.render.render_debra_png.

Rendering is exercised on a small synthetic uint8 composite over the shared
west-Texas lat/lon fixture; assertions are structural (a real PNG of sensible
size lands where asked), not pixel-exact. Cartopy's Natural Earth cache is
already populated by the quicklook path, so no network is needed.
"""

import datetime as dt

import matplotlib.image as mpimg
import numpy as np
import pytest

from shachen.render import render_debra_png

WHEN = dt.datetime(2017, 3, 23, 21, 45)


def _rgb(ny: int = 20, nx: int = 20) -> np.ndarray:
    ramp = np.linspace(0, 255, ny * nx, dtype=np.float64).reshape(ny, nx)
    rgb = np.stack([ramp, ramp[::-1], np.full((ny, nx), 32.0)], axis=-1)
    return rgb.astype(np.uint8)


def test_writes_png(tmp_path, latlon_area):
    out = tmp_path / "debra.png"
    result = render_debra_png(_rgb(), latlon_area, WHEN, out)
    assert result == out
    assert out.exists()
    data = out.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    img = mpimg.imread(out)
    assert img.ndim == 3
    assert img.shape[0] >= 100 and img.shape[1] >= 100


def test_custom_title_and_dpi(tmp_path, latlon_area):
    out = tmp_path / "titled.png"
    result = render_debra_png(_rgb(), latlon_area, WHEN, out, title="demo", dpi=72)
    assert result.exists()
    assert result.stat().st_size > 0


@pytest.mark.parametrize(
    "bad",
    [
        _rgb().astype(np.float64),  # wrong dtype
        _rgb()[..., 0],  # missing channel axis
        np.zeros((20, 20, 4), dtype=np.uint8),  # RGBA not RGB
    ],
)
def test_rejects_bad_rgb(tmp_path, latlon_area, bad):
    with pytest.raises(ValueError):
        render_debra_png(bad, latlon_area, WHEN, tmp_path / "bad.png")
