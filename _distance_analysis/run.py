"""
run.py
------
Entry point: compute centroid → bottom-edge distances and save CSVs.

Edit ``nuclear_distance_analysis/config.py`` to set your input/output
paths before running.

Usage
-----
    python run.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from skimage import io

from nuclear_distance_analysis import (
    analyze_new_mask,
    create_beeswarm_boxplot,
    find_merged_masks_in_subfolders,
    key_from_bottom,
    key_from_mMask,
    key_from_original,
    list_images,
    list_originals_strict,
    load_as_2d,
    long_df_from_groups,
    make_new_mask,
    print_file_summary,
    safe_slug,
)
from nuclear_distance_analysis.config import (
    BOTTOM_DIR,
    MERGED_DIR,
    MMASK_DIR,
    ORIG_DIR,
    ORIG_SUFFIX,
    OUT_DIR,
)


# =============================================================================
# Output directories
# =============================================================================
OUT_NEW_MASK_DIR = OUT_DIR / "new_masks_merged_inside_mMask"
OUT_CSV_DIR      = OUT_DIR / "csv"
OUT_PLOTS_DIR    = OUT_DIR / "plots"

for d in (OUT_DIR, OUT_NEW_MASK_DIR, OUT_CSV_DIR, OUT_PLOTS_DIR):
    d.mkdir(parents=True, exist_ok=True)


# =============================================================================
# File discovery and matching
# =============================================================================
orig_files   = list_originals_strict(ORIG_DIR, ORIG_SUFFIX)
mMask_files   = list_images(MMASK_DIR)
bottom_files = list_images(BOTTOM_DIR)
merged_map   = find_merged_masks_in_subfolders(MERGED_DIR)

orig_map   = {key_from_original(p.name): p for p in orig_files}
mMask_map   = {key_from_mMask(p.name):     p for p in mMask_files}
bottom_map = {key_from_bottom(p.name):   p for p in bottom_files}

all_keys = sorted(set(orig_map) & set(merged_map) & set(mMask_map) & set(bottom_map))
print_file_summary(orig_map, merged_map, mMask_map, bottom_map, all_keys)


# =============================================================================
# Batch loop
# =============================================================================
all_components, all_label_summaries, sample_counts_rows = [], [], []

for k in all_keys:
    orig   = load_as_2d(orig_map[k])
    merged = load_as_2d(merged_map[k]).astype(np.int32)
    mMask   = load_as_2d(mMask_map[k]).astype(np.uint8)
    bottom = load_as_2d(bottom_map[k]).astype(np.uint8)

    new_mask = make_new_mask(merged, mMask)

    io.imsave(
        str(OUT_NEW_MASK_DIR / f"{safe_slug(k)}_newmask_merged_inside_mMask.tif"),
        new_mask.astype(np.int32),
        check_contrast=False,
    )

    df_comp, df_lab, _y_edge, sample_counts = analyze_new_mask(new_mask, bottom, k)

    for df_part in (df_comp, df_lab):
        if not df_part.empty:
            df_part["object_count_components_sample"] = sample_counts["object_count_components"]
            df_part["object_count_labels_sample"]     = sample_counts["object_count_labels"]

    df_comp.to_csv(OUT_CSV_DIR / f"{safe_slug(k)}_components_distances.csv",    index=False)
    df_lab.to_csv( OUT_CSV_DIR / f"{safe_slug(k)}_label_summary_distances.csv", index=False)

    all_components.append(df_comp)
    all_label_summaries.append(df_lab)
    sample_counts_rows.append(sample_counts)


# =============================================================================
# Global CSVs
# =============================================================================
df_all_components      = pd.concat(all_components,      ignore_index=True) if all_components      else pd.DataFrame()
df_all_label_summaries = pd.concat(all_label_summaries, ignore_index=True) if all_label_summaries else pd.DataFrame()
df_sample_counts       = pd.DataFrame(sample_counts_rows)

df_all_components.to_csv(     OUT_CSV_DIR / "ALL_components_distances.csv",    index=False)
df_all_label_summaries.to_csv(OUT_CSV_DIR / "ALL_label_summary_distances.csv", index=False)
df_sample_counts.to_csv(      OUT_CSV_DIR / "ALL_sample_counts.csv",           index=False)


# =============================================================================
# Per-sample mean of dist_mean_px  (distance to bottom edge — MEAN)
# =============================================================================
if not df_all_label_summaries.empty:
    df_mean_summary = (
        df_all_label_summaries
        .groupby("sample")["dist_mean_px"]
        .mean()
        .reset_index()
        .rename(columns={"dist_mean_px": "dist_to_bottom_edge_mean_px"})
    )
    df_mean_summary.to_csv(OUT_CSV_DIR / "distance_to_bottom_edge_mean.csv", index=False)

    print("\nDistance to Bottom Edge — MEAN (per sample):")
    print(df_mean_summary.to_string(index=False))

    # Beeswarm + boxplot — MEAN
    samples = sorted(df_all_label_summaries["sample"].dropna().unique())
    groups_mean = {
        s: df_all_label_summaries.loc[df_all_label_summaries["sample"] == s, "dist_mean_px"]
           .dropna().to_numpy(float)
        for s in samples
    }
    df_mean_plot = long_df_from_groups(groups_mean, y_name="dist_mean_px")

    counts_map = (
        df_sample_counts.set_index("sample")[["object_count_components", "object_count_labels"]]
        .to_dict("index")
        if not df_sample_counts.empty else {}
    )
    tick_labels = [
        f"{s}\nNcomp={counts_map.get(s, {}).get('object_count_components', 0)}, "
        f"Nlab={counts_map.get(s, {}).get('object_count_labels', 0)}"
        for s in samples
    ]

    create_beeswarm_boxplot(
        df_mean_plot, y_column="dist_mean_px", order=samples, tick_labels=tick_labels,
        title="Distance to Bottom Edge (Beeswarm) — MEAN",
        swarm_color="mediumblue", figsize=(max(12, len(samples) * 2.2), 6),
        save_png=OUT_PLOTS_DIR / "beeswarm_boxplot_mean.png",
    )


print("\n✅ Processing complete.")
print(f"   Components     : {df_all_components.shape}")
print(f"   Label summaries: {df_all_label_summaries.shape}")
print(f"   Sample counts  : {df_sample_counts.shape}")
print(f"   CSVs saved to  : {OUT_CSV_DIR}")
