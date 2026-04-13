from typing import Tuple
"""
analysis.py
-----------
Per-label / per-component centroid-to-edge distance measurements.
"""

import numpy as np
import pandas as pd
from skimage.measure import label as cc_label, regionprops

from .masks import compute_bottom_edge_profile, distance_centroid_to_bottom_edge


def analyze_new_mask(
    new_mask: np.ndarray,
    bottom_mask: np.ndarray,
    sample_key: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, dict]:
    """
    Compute centroid-to-edge distances for every connected component in *new_mask*.

    Parameters
    ----------
    new_mask    : labelled mask (merged ∩ mMask)
    bottom_mask : binary bottom-tissue mask used to derive the edge profile
    sample_key  : string identifier for this sample

    Returns
    -------
    df_components    : one row per connected component
    df_label_summary : one row per label — aggregated distance statistics
    y_edge           : bottom-edge profile (pixels, per column)
    sample_counts    : dict with ``object_count_labels`` and ``object_count_components``
    """
    y_edge = compute_bottom_edge_profile(bottom_mask)
    labels = np.unique(new_mask)
    labels = labels[labels > 0]

    comp_rows, summary_rows = [], []
    total_components = 0

    for lab in labels:
        cc    = cc_label(new_mask == lab, connectivity=1)
        props = regionprops(cc)
        if not props:
            continue

        dists = []
        for i, rp in enumerate(props, start=1):
            dist = distance_centroid_to_bottom_edge(rp.centroid, y_edge)
            dists.append(dist)
            comp_rows.append({
                "sample":                    sample_key,
                "label_id":                  int(lab),
                "component_id_within_label": int(i),
                "area_pixels":               int(rp.area),
                "centroid_y":                float(rp.centroid[0]),
                "centroid_x":                float(rp.centroid[1]),
                "dist_to_edge_bottom_px":    float(dist),
            })

        total_components += len(dists)
        d = np.asarray(dists, dtype=float)
        summary_rows.append({
            "sample":         sample_key,
            "label_id":       int(lab),
            "n_components":   int(len(dists)),
            "dist_mean_px":   float(np.nanmean(d)),
            "dist_median_px": float(np.nanmedian(d)),
            "dist_min_px":    float(np.nanmin(d)),
            "dist_max_px":    float(np.nanmax(d)),
        })

    sample_counts = {
        "sample":                  sample_key,
        "object_count_labels":     int(labels.size),
        "object_count_components": int(total_components),
    }
    return (
        pd.DataFrame(comp_rows),
        pd.DataFrame(summary_rows),
        y_edge,
        sample_counts,
    )
