"""
masks.py
--------
Mask construction and bottom-edge geometry utilities.
"""

import numpy as np


def make_new_mask(merged_mask: np.ndarray, mMask_mask: np.ndarray) -> np.ndarray:
    """Zero out everything in *merged_mask* that lies outside the mMask region."""
    new_mask = merged_mask.copy()
    new_mask[mMask_mask <= 0] = 0
    return new_mask


def compute_bottom_edge_profile(bottom_mask_2d: np.ndarray) -> np.ndarray:
    """
    Return the y-coordinate of the lowest foreground pixel in each column.
    Columns with no foreground are filled by linear interpolation.
    """
    m = bottom_mask_2d > 0
    _, w = m.shape
    y_edge = np.full(w, np.nan)

    for x in range(w):
        ys = np.where(m[:, x])[0]
        if ys.size:
            y_edge[x] = ys.max()

    valid = np.where(~np.isnan(y_edge))[0]
    if valid.size:
        y_edge = np.interp(np.arange(w), valid, y_edge[valid])
    return y_edge


def distance_centroid_to_bottom_edge(
    centroid_yx: tuple,
    y_edge: np.ndarray,
) -> float:
    """Signed pixel distance from a centroid to the bottom edge at the same column."""
    cy, cx = centroid_yx
    if y_edge is None or np.all(np.isnan(y_edge)):
        return np.nan
    x0 = int(np.clip(round(cx), 0, len(y_edge) - 1))
    return float(y_edge[x0] - cy)
