"""Shared fixtures: small pyresample areas for the acceptance tests."""

import pytest
from pyresample.geometry import AreaDefinition


def _latlon_area(
    lon_min: float, lat_min: float, lon_max: float, lat_max: float, width: int, height: int
) -> AreaDefinition:
    return AreaDefinition(
        "test_latlon",
        "regular lat/lon test grid",
        "test_latlon",
        {"proj": "longlat", "datum": "WGS84"},
        width,
        height,
        (lon_min, lat_min, lon_max, lat_max),
    )


@pytest.fixture
def make_latlon_area():
    """Factory for small regular lat/lon AreaDefinitions."""
    return _latlon_area


@pytest.fixture
def latlon_area() -> AreaDefinition:
    """20x20 grid over 30-34N, 104-100W (west Texas / SE New Mexico, all land)."""
    return _latlon_area(-104.0, 30.0, -100.0, 34.0, 20, 20)


@pytest.fixture
def geos_area() -> AreaDefinition:
    """GOES-16-like full disk at 50x50 pixels; corner pixels fall off the disk."""
    return AreaDefinition(
        "test_geos",
        "GOES-16-like full disk",
        "test_geos",
        {
            "proj": "geos",
            "h": 35786023.0,
            "lon_0": -89.5,
            "sweep": "x",
            "a": 6378137.0,
            "b": 6356752.31414,
        },
        50,
        50,
        (-5434894.885, -5434894.885, 5434894.885, 5434894.885),
    )
