"""The Eq. 3 normalization primitive used by every DEBRA test."""

import numpy as np

from shachen.constants import Bounds


def normalize(x, bounds: Bounds):
    """N(x) = clip((x - MIN) / (MAX - MIN), 0, 1)  (Miller et al. 2017, Eq. 3).

    Works on scalars, numpy arrays, and xarray DataArrays. ``bounds.min`` may
    exceed ``bounds.max`` (used for the cos-zenith blends specified in
    zenith-angle space): the sense of the ramp reverses.
    """
    return np.clip((x - bounds.min) / (bounds.max - bounds.min), 0.0, 1.0)


def normalize_cos_zenith(zenith_deg, zenith_bounds_deg: Bounds, exponent: float):
    """N(cos(theta); cos(MIN_deg), cos(MAX_deg)) ** exponent (Eqs. 20-21, 25)."""
    cos_bounds = Bounds(
        np.cos(np.deg2rad(zenith_bounds_deg.min)),
        np.cos(np.deg2rad(zenith_bounds_deg.max)),
    )
    return normalize(np.cos(np.deg2rad(zenith_deg)), cos_bounds) ** exponent
