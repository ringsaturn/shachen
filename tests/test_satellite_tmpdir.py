"""``load_scene`` must not leave satpy's unpacked bz2 segments behind.

satpy's AHI HSD handler unpacks each ``.bz2`` segment into ``config["tmp_dir"]``
and registers ``weakref.finalize(self, self._cleanup)`` to delete it -- a bound
method that keeps the handler alive, so the file survives the Scene and is
only removed at interpreter exit (never, when a worker is killed). Measured on
2026-09-17: 93 GB of leaked segments on one batch worker. ``load_scene``
therefore gives every call its own temporary directory and computes the
arrays before it is removed. This test stands in a fake ``Scene`` that
behaves like the leaky handler: it writes a file into satpy's tmp_dir during
``load`` and never deletes it.
"""

from pathlib import Path

import numpy as np
import pytest
import satpy
import xarray as xr
from pyresample import create_area_def

import shachen.io.satellite as satellite
from shachen.constants import Band


class _LeakyScene:
    """The parts of ``satpy.Scene`` that ``load_scene`` touches, leaking a file."""

    leaked: list[Path] = []

    def __init__(self, filenames, reader):
        self._area = create_area_def(
            "t", {"proj": "latlong"}, shape=(4, 6), area_extent=(0, 0, 6, 4)
        )
        self._data = {}

    def load(self, names):
        tmp_dir = Path(satpy.config.get("tmp_dir"))
        leaked = tmp_dir / "segment.DAT"
        leaked.write_bytes(b"x" * 16)
        type(self).leaked.append(leaked)
        for name in names:
            self._data[name] = xr.DataArray(
                np.full((4, 6), 280.0), dims=("y", "x"), attrs={"area": self._area}
            ).chunk()

    def coarsest_area(self):
        return self._area

    def resample(self, area, resampler):
        return self

    def crop(self, ll_bbox):
        return self

    @property
    def start_time(self):
        return "2026-09-17T00:00:00"

    def __getitem__(self, name):
        return self._data[name]


def test_load_scene_leaves_no_temp_files(monkeypatch, tmp_path):
    monkeypatch.setattr(satellite, "Scene", _LeakyScene)
    _LeakyScene.leaked.clear()
    # Whatever satpy's tmp_dir was before, load_scene must not write there.
    with satpy.config.set(tmp_dir=str(tmp_path)):
        ds = satellite.load_scene(["a.bz2"], reader="ahi_hsd", roles=[Band.TIR_104])

    assert _LeakyScene.leaked, "the fake scene did not exercise tmp_dir"
    assert not any(p.exists() for p in _LeakyScene.leaked)
    assert "shachen-l1b-" in str(_LeakyScene.leaked[0])
    assert not list(tmp_path.iterdir())
    # The arrays were computed inside the sandbox: numpy, not a dask graph
    # that would read the (now deleted) segment on first use.
    assert isinstance(ds["bt_tir_104"].data, np.ndarray)
    assert float(ds["bt_tir_104"].mean()) == pytest.approx(280.0)
