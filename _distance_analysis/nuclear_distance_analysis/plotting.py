"""
plotting.py
-----------
Beeswarm + boxplot figure for the distance-to-bottom-edge mean.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — no display required
import matplotlib.pyplot as plt
import seaborn as sns


def long_df_from_groups(
    groups_dict: Dict[str, np.ndarray],
    y_name: str,
    group_name: str = "sample",
) -> pd.DataFrame:
    """Convert ``{sample: values_array}`` to a long-format DataFrame."""
    rows = []
    for s, vals in groups_dict.items():
        vals = np.asarray(vals, dtype=float)
        for v in vals[~np.isnan(vals)]:
            rows.append({group_name: s, y_name: float(v)})
    return pd.DataFrame(rows)


def create_beeswarm_boxplot(
    df: pd.DataFrame,
    y_column: str,
    group_column: str = "sample",
    title: str = "Distance to Bottom Edge",
    figsize: tuple = (18, 6),
    swarm_color: str = "darkgreen",
    swarm_alpha: float = 0.5,
    order: Optional[List] = None,
    tick_labels: Optional[List] = None,
    decimals: int = 5,
    headroom_frac: float = 0.18,
    text_base_frac: float = 0.02,
    text_gap_frac: float = 0.05,
    save_png: Optional[Path] = None,
    save_pdf: Optional[Path] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Beeswarm + boxplot (white boxes, black outlines, no outlier markers).
    Annotates per-group mean and median above each box.
    Saves to disk; never displays interactively.
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.set_style("whitegrid")

    if order is None:
        order = sorted(df[group_column].dropna().unique())

    sns.boxplot(
        data=df, x=group_column, y=y_column, order=order, ax=ax,
        color="white", width=0.55, linewidth=2.0, fliersize=0, showmeans=False,
        boxprops=dict(edgecolor="black", linewidth=2.0),
        whiskerprops=dict(color="black", linewidth=2.0),
        capprops=dict(color="black", linewidth=2.0),
        medianprops=dict(color="black", linewidth=2.5),
    )
    sns.swarmplot(
        data=df, x=group_column, y=y_column, order=order, ax=ax,
        color=swarm_color, alpha=swarm_alpha, size=3, edgecolor="none",
    )

    ax.set_xlabel("Sample (Batch)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Distance (pixels)", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")

    if tick_labels is not None:
        ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    else:
        plt.xticks(rotation=45, ha="right")

    ax.yaxis.grid(True, alpha=0.3, linestyle="--", zorder=0)
    ax.set_axisbelow(True)

    # Annotate mean and median above each box
    y_min, y_max = ax.get_ylim()
    y_range = (y_max - y_min) or 1.0
    ax.set_ylim(y_min, y_max + headroom_frac * y_range)
    y_base = y_max + text_base_frac * y_range
    y_gap  = text_gap_frac * y_range

    for i, grp in enumerate(order):
        vals = df.loc[df[group_column] == grp, y_column].dropna()
        if vals.empty:
            continue
        ax.text(i, y_base,         f"Med={vals.median():.{decimals}f}",
                ha="center", va="bottom", fontsize=10, fontweight="bold",
                color="black", clip_on=False)
        ax.text(i, y_base + y_gap, f"Mean={vals.mean():.{decimals}f}",
                ha="center", va="bottom", fontsize=9, style="italic",
                color="darkblue", clip_on=False)

    plt.tight_layout()
    if save_png:
        fig.savefig(save_png, dpi=300, bbox_inches="tight")
    if save_pdf:
        fig.savefig(save_pdf, bbox_inches="tight")
    plt.close(fig)
    return fig, ax
