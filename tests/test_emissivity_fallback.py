"""CAMEL month fallback: what happens after the product stops being made."""

import datetime as dt
from pathlib import Path

import pytest

from shachen.io import emissivity as emis


@pytest.fixture
def archive(tmp_path, monkeypatch):
    """A fake CAMEL archive: only the months listed exist, and a 'download'
    just writes the granule name into the scratch file."""

    def make(months: set[str]):
        def search(first):
            return [f"CAM5K30EM_emis_{first:%Y%m}_V003.nc"] if f"{first:%Y%m}" in months else []

        def download(matched, destination):
            path = Path(destination) / matched[0]
            path.write_text(matched[0])
            return [str(path)]

        monkeypatch.setattr(emis, "_search_month", search)
        fake = type(
            "earthaccess",
            (),
            {"login": staticmethod(lambda **kw: None), "download": staticmethod(download)},
        )
        monkeypatch.setitem(__import__("sys").modules, "earthaccess", fake)
        return tmp_path

    return make


def test_the_month_itself_wins_when_it_exists(archive):
    out = archive({"202304"})
    path = emis.fetch_emissivity(dt.date(2023, 4, 15), out)
    assert path.name == "CAM5K30EM_202304.nc"


def test_a_month_past_the_products_end_falls_back_a_year(archive, recwarn):
    """CAM5K30EM V003 stops at 2023-12; an April 2025 scene still needs April."""
    out = archive({"202304"})
    path = emis.fetch_emissivity(dt.date(2025, 4, 11), out)
    assert path.name == "CAM5K30EM_202304.nc"
    assert "no 2025-04 granule" in str(recwarn.pop(UserWarning).message)


def test_the_fallback_never_crosses_into_another_month(archive):
    """March emissivity is not April's: only whole years are walked back."""
    out = archive({"202303"})
    with pytest.raises(RuntimeError, match="same month of the"):
        emis.fetch_emissivity(dt.date(2025, 4, 11), out, fallback_years=3)


def test_a_cached_fallback_is_reused_without_searching(archive, recwarn):
    out = archive({"202304"})
    emis.fetch_emissivity(dt.date(2025, 4, 11), out)
    recwarn.clear()
    # Now with an empty archive: only the cache can answer.
    out = archive(set())
    path = emis.fetch_emissivity(dt.date(2025, 4, 11), out)
    assert path.name == "CAM5K30EM_202304.nc"
    assert "using 2023-04" in str(recwarn.pop(UserWarning).message)
