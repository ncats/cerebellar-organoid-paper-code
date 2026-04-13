# nuclear_distance_analysis

Measures the distance of segmented nuclei centroids to the bottom tissue edge and saves results as CSV and png files.

## Setup

```bash
pip install -r requirements.txt
```

## Configuration

Edit `nuclear_distance_analysis/config.py` and set the five path variables:

| Variable | Description |
|---|---|
| `ORIG_DIR` | Folder containing original `.tif` images |
| `MERGED_DIR` | Root folder with one subfolder per sample (merged masks) |
| `MMASK_DIR` | Folder containing mMask files |
| `BOTTOM_DIR` | Folder containing bottom mask files |
| `OUT_DIR` | Output folder for CSVs and new masks |

## Usage

```bash
python run.py
```

## Output

All files are written to `OUT_DIR/`:

```
OUT_DIR/
├── csv/
│   ├── ALL_components_distances.csv       # one row per connected component
│   ├── ALL_label_summary_distances.csv    # one row per label (mean, median, min, max)
│   ├── ALL_sample_counts.csv              # object counts per sample
│   ├── distance_to_bottom_edge_mean.csv   # per-sample mean distance (beeswarm — MEAN)
│   └── <sample>_*.csv                     # per-sample CSVs
└── new_masks_merged_inside_mMask/
    └── <sample>_newmask_merged_inside_mMask.tif
```

### Key output: `distance_to_bottom_edge_mean.csv`

| Column | Description |
|---|---|
| `sample` | Sample identifier |
| `dist_to_bottom_edge_mean_px` | Mean of per-label mean distances (pixels) |
