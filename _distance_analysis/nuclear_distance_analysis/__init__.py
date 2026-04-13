"""
nuclear_distance_analysis
=========================
Pipeline: distance of segmented nuclei to the bottom tissue edge.
"""

from .analysis import analyze_new_mask
from .io_utils import (
    find_merged_masks_in_subfolders,
    key_from_bottom,
    key_from_mMask,
    key_from_original,
    list_images,
    list_originals_strict,
    load_as_2d,
    print_file_summary,
    safe_slug,
)
from .masks import (
    compute_bottom_edge_profile,
    distance_centroid_to_bottom_edge,
    make_new_mask,
)
from .plotting import (
    create_beeswarm_boxplot,
    long_df_from_groups,
)

__all__ = [
    # analysis
    "analyze_new_mask",
    # io
    "find_merged_masks_in_subfolders",
    "key_from_bottom",
    "key_from_mMask",
    "key_from_original",
    "list_images",
    "list_originals_strict",
    "load_as_2d",
    "print_file_summary",
    "safe_slug",
    # masks
    "compute_bottom_edge_profile",
    "distance_centroid_to_bottom_edge",
    "make_new_mask",
    # plotting
    "create_beeswarm_boxplot",
    "long_df_from_groups",
]
