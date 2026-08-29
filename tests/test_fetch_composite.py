"""Acceptance tests for the composite-window additions to scripts/fetch_case.py
: the pure day enumeration and the channel/constant wiring.
The S3 download itself is exercised manually (network)."""

import datetime as dt
import inspect
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import fetch_case  # noqa: E402

from shachen.composite import COMPOSITE_BANDS  # noqa: E402
from shachen.constants import (  # noqa: E402
    ABI_BANDS,
    COMPOSITE_MIN_DAYS,
    COMPOSITE_WINDOW_DAYS,
)

_WHEN = dt.datetime(2020, 12, 13, 21, 20)


def test_window_constants():
    assert COMPOSITE_WINDOW_DAYS == 14
    assert COMPOSITE_MIN_DAYS == 5


def test_composite_channels_match_composite_bands():
    assert fetch_case.COMPOSITE_CHANNELS == ("C11", "C13", "C15")
    assert fetch_case.COMPOSITE_CHANNELS == tuple(ABI_BANDS[band] for band in COMPOSITE_BANDS)


def test_composite_days_precede_case_day_chronologically():
    days = fetch_case.composite_days(_WHEN, 14)
    assert len(days) == 14
    assert days[0] == dt.date(2020, 11, 29)  # crosses the month boundary
    assert days[-1] == dt.date(2020, 12, 12)
    assert all(a < b for a, b in zip(days, days[1:], strict=False))
    assert _WHEN.date() not in days


def test_composite_days_single():
    assert fetch_case.composite_days(_WHEN, 1) == [dt.date(2020, 12, 12)]


def test_fetch_composite_default_window():
    default = inspect.signature(fetch_case.fetch_composite).parameters["days"].default
    assert default == COMPOSITE_WINDOW_DAYS
