"""Unit tests for the Eq. 3 normalization primitive."""

import numpy as np
import pytest

from shachen.constants import Bounds
from shachen.norm import normalize, normalize_cos_zenith, normalize_interp


def test_scalar_ramp():
    b = Bounds(0.0, 10.0)
    assert normalize(-5.0, b) == 0.0
    assert normalize(0.0, b) == 0.0
    assert normalize(5.0, b) == 0.5
    assert normalize(10.0, b) == 1.0
    assert normalize(25.0, b) == 1.0


def test_array_input():
    b = Bounds(2.0, 4.5)  # the CM3 bounds
    x = np.array([0.0, 2.0, 3.25, 4.5, 10.0])
    np.testing.assert_allclose(normalize(x, b), [0.0, 0.0, 0.5, 1.0, 1.0])


def test_reversed_bounds_reverse_the_ramp():
    b = Bounds(10.0, 0.0)
    assert normalize(0.0, b) == 1.0
    assert normalize(10.0, b) == 0.0
    assert normalize(5.0, b) == 0.5


def test_cos_zenith_blend_endpoints():
    # Eq. 21 (erratum): B_trm,day ramps over zenith 90 deg -> 75 deg
    b = Bounds(90.0, 75.0)
    assert normalize_cos_zenith(95.0, b, 1.5) == 0.0  # night side: fully off
    assert normalize_cos_zenith(90.0, b, 1.5) == 0.0
    assert normalize_cos_zenith(75.0, b, 1.5) == pytest.approx(1.0)
    assert normalize_cos_zenith(30.0, b, 1.5) == 1.0  # full day: fully on
    mid = normalize_cos_zenith(82.0, b, 1.5)
    assert 0.0 < mid < 1.0


def test_cos_zenith_exponent_applied():
    b = Bounds(90.0, 75.0)
    lin = normalize_cos_zenith(82.0, b, 1.0)
    assert normalize_cos_zenith(82.0, b, 1.5) == pytest.approx(lin**1.5)


def test_interp_endpoints_pick_one_interval_or_the_other():
    day = Bounds(0.40, 2.50)
    ngt = Bounds(0.20, 1.25)
    assert normalize_interp(1.45, day, ngt, 1.0) == pytest.approx(normalize(1.45, day))
    assert normalize_interp(1.45, day, ngt, 0.0) == pytest.approx(normalize(1.45, ngt))


def test_interp_midpoint_uses_the_averaged_interval():
    day = Bounds(0.40, 2.50)
    ngt = Bounds(0.20, 1.25)
    mid = Bounds((0.40 + 0.20) / 2, (2.50 + 1.25) / 2)
    assert normalize_interp(1.0, day, ngt, 0.5) == pytest.approx(normalize(1.0, mid))


def test_interp_accepts_per_element_weights():
    day = Bounds(0.0, 2.0)
    ngt = Bounds(0.0, 1.0)
    weight = np.array([0.0, 1.0])
    np.testing.assert_allclose(
        normalize_interp(np.array([0.5, 0.5]), day, ngt, weight), [0.5, 0.25]
    )


def test_interp_clips_to_the_unit_interval():
    day = Bounds(0.40, 2.50)
    ngt = Bounds(0.20, 1.25)
    assert normalize_interp(-3.0, day, ngt, 0.3) == 0.0
    assert normalize_interp(99.0, day, ngt, 0.3) == 1.0
