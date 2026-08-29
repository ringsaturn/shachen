"""Unit tests for the Eq. 3 normalization primitive."""

import numpy as np
import pytest

from shachen.constants import Bounds
from shachen.norm import normalize, normalize_cos_zenith


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
