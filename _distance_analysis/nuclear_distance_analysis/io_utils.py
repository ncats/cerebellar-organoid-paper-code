from typing import Dict, List
"""
io_utils.py
-----------
File discovery, sample-key extraction, and image loading utilities.
"""

import re
from pathlib import Path

import numpy as np
from skimage import io


# =============================================================================
# Helpers
# =============================================================================

def safe_slug(s: str) -> str:
    """Replace non-alphanumeric characters with underscores."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", s)


def list_images(folder: Path, exts: tuple = (".tif", ".tiff")) -> List[Path]:
    """Return all image files in *folder* matching *exts* (case-insensitive)."""
    files = []
    for ext in exts:
        files.extend(folder.glob(f"*{ext}"))
        files.extend(folder.glob(f"*{ext.upper()}"))
    return sorted(set(files))


def list_originals_strict(folder: Path, required_suffix: str) -> List[Path]:
    """Return only images whose name ends with *required_suffix* (case-insensitive)."""
    suf     = required_suffix.lower()
    suf_alt = (suf[:-4] + ".tiff") if suf.endswith(".tif") else None

    picked = []
    for p in list_images(folder):
        n = p.name.lower()
        if n.endswith(suf) or (suf_alt and n.endswith(suf_alt)):
            picked.append(p)
    return sorted(picked)


# =============================================================================
# Sample-key extractors
# =============================================================================

def key_from_original(name: str) -> str:
    return name.split("___")[0].strip() if "___" in name else Path(name).stem.strip()


def key_from_mMask(name: str) -> str:
    if "___" in name:
        return name.split("___")[0].strip()
    return Path(name).stem.replace("_mMask_mask", "").strip()


def key_from_bottom(name: str) -> str:
    if "___" in name:
        return name.split("___")[0].strip()
    return Path(name).stem.replace("_bottom_mask", "").strip()


# =============================================================================
# Merged-mask discovery
# =============================================================================

def find_merged_masks_in_subfolders(merged_root: Path) -> Dict[str, Path]:
    """
    Expect::

        merged_root/
            <sample_key>/
                <sample_key>_merged_objects.tif[f]

    Returns ``{sample_key: Path}`` for every matching subfolder.
    """
    merged_map: Dict[str, Path] = {}
    for sub in merged_root.iterdir():
        if not sub.is_dir():
            continue
        key = sub.name.strip()
        for suffix in ("_merged_objects.tif", "_merged_objects.tiff"):
            candidate = sub / f"{key}{suffix}"
            if candidate.exists():
                merged_map[key] = candidate
                break
        else:
            candidates = (
                list(sub.glob("*_merged_objects.tif")) +
                list(sub.glob("*_merged_objects.tiff"))
            )
            if candidates:
                exact = [c for c in candidates if c.name.startswith(key)]
                merged_map[key] = exact[0] if exact else candidates[0]
    return merged_map


# =============================================================================
# Image loading
# =============================================================================

def load_as_2d(path: Path) -> np.ndarray:
    """Load any 2-D / 3-D TIFF and return a single 2-D plane."""
    arr = io.imread(str(path))
    if arr.ndim == 2:
        return arr
    if arr.ndim == 3:
        return arr[..., 0] if arr.shape[-1] <= 4 else arr[0]
    while arr.ndim > 2:
        arr = arr[0]
    return arr


# =============================================================================
# Reporting
# =============================================================================

def print_file_summary(
    orig_map: dict,
    merged_map: dict,
    mMask_map: dict,
    bottom_map: dict,
    all_keys: list,
) -> None:
    sections = [
        ("Original files", list(orig_map.values())),
        ("Merged masks",   list(merged_map.values())),
        ("MMask masks",     list(mMask_map.values())),
        ("Bottom masks",   list(bottom_map.values())),
    ]
    for title, files in sections:
        print(f"\n{title} ({len(files)} files)")
        print("-" * (len(title) + 12))
        for f in files:
            print(f)

    print(f"\nMatched samples : {len(all_keys)}")
    if not all_keys:
        print("⚠️  No matches found. Sample key previews:")
        for label, d in (("orig", orig_map), ("merged", merged_map),
                          ("mMask", mMask_map), ("bottom", bottom_map)):
            print(f"  {label:8s}: {list(d)[:5]}")
