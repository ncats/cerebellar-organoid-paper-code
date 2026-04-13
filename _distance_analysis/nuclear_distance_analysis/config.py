"""
config.py
---------
All user-facing settings in one place.
Edit the five path variables below before running.
"""

from pathlib import Path

# =============================================================================
# REQUIRED: set these five paths before running
# =============================================================================

ORIG_DIR_path = r"./data/ORIG_IMGS"
MERGED_DIR_path = r"./data/MERGED_"
MMASK_DIR_path = r"./data/MASK_"
BOTTOM_DIR_path = r"./data/BOTTOM_"
TUJ1_QUANTIFICATION_path = r"./data/TUJ1_QUANTIFICATION_"

ORIG_DIR   = Path(ORIG_DIR_path)
MERGED_DIR = Path(MERGED_DIR_path)
MMASK_DIR   = Path(MMASK_DIR_path)
BOTTOM_DIR = Path(BOTTOM_DIR_path)
OUT_DIR    = Path(TUJ1_QUANTIFICATION_path)

# Only originals whose filename ends with this suffix are processed.
ORIG_SUFFIX = "___Ch1-T2__DAPI.tif"
