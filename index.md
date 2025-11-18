# read\_image\_patch\_colors
# Reference and Usage Guide
**Version:** 1.4<br>
**Author:** Knut Larsson<br>
**Purpose:** Generate ArgyllCMS compatible `.ti1` and `.ti2` files from a grid of color patches in an image.

## Table of Contents
- [Overview](#overview)
- [Image Input](#image-input)
- [Grid Geometry](#grid-geometry)
- [Sampling](#sampling)
- [Labeling Rules](#labeling-rules)
- [Color Processing](#color-processing)
- [Output Files](#output-files)
- [Patch Value Ranges](#patch-value-ranges)
- [Valid Color Blocks](#valid-color-blocks)
- [Dependencies](#dependencies)
- [Command-line Arguments](#command-line-arguments)
- [Notes on Accuracy](#notes-on-accuracy)
- [Patch Sampling Method](#patch-sampling-method)
- [Output File Format (TI2)](#output-file-format-ti2)
- [Error Conditions](#error-conditions)
- [Workflow Summary](#workflow-summary)
- [Examples](#examples)

# Overview
This program extracts color values from a rectangular grid of color patches in an image, computes colorimetric values (RGB percentages, XYZ, Lab), applies row/column labeling rules, and writes three output files:

- ArgyllCMS .ti1 file
- ArgyllCMS .ti2 file
- CSV file (space-separated)

It supports sampling in mean or median mode, configurable patch geometry, numeric or alphabetic labeling patterns, and RGB/XYZ/Lab output combinations.

The script is designed for use in color-management workflows where printed color targets must be scanned or photographed and converted into Argyll measurement files.

Example use cases: 

1. Using the image of a reference target, create a `.ti1` file so thatn one may use ArgyllCMS `printtarg` command and generate a target using the colors of the image, which then can be used to create a printer profile.

2. Using the image of a reference target, create a `.ti2` file so that one may print the target, scan it, and then use ArgyllCMS `scanin` command to create a printer profile.

# Image Input
Accepted image types: any PIL-compatible raster file
Color depth: handled as 16-bit RGB internally
Color space:

- "srgb" (default)
- "adobergb"

Parameter: `--image` / `-i`
Provides the path to the patch-grid image.

# Grid Geometry
Geometry is computed from:

- `--patch_first_xy` — top-left patch center (x,y)
- `--patch_last_xy` — bottom-right patch center (x,y)
- `--num_cols` — number of columns
- `--num_rows` — number of rows
- `--patch_width_height_ratio` — width/height of a patch

The script calculates patch center coordinates and boundaries.

# Sampling
Parameter: `--sample_fraction`
Range: `0 < f ≤ 0.6`
Defines fraction of the *smaller* patch dimension used for sampling.

Sampling modes:

- “mean”
- “median”

Sampling region is square, centered in each patch.

# Labeling Rules
Row and column labels come from:

- `--row_labels` (e.g., "A-D", "1-20")
- `--col_labels`

Accepted patterns:

- Alphabetic: A, B, … Z, AA, AB, …
- Numeric: integers ≥ 0

Label order:

- `--patch_label_order = col_then_row`
- `--patch_label_order = row_then_col`

Label counts must match grid dimensions exactly.

# Color Processing
Sampled RGB → RGB percentage (0–100):
`RGB_percent = (R_16bit / 65535) * 100`

Conversions use python-colormath according to selected color space (“srgb” / “adobergb”).

The patch with the highest XYZ_Y is stored as `APPROX_WHITE_POINT`.

XYZ output:
• CIE 1931 2°
• Chromatically adapted to D50
• White scaled to Y=100

Lab output:
• CIE Lab (1976)
• D50 reference

# Output Files
Three files:

### 1. TI1 (.ti1)
Contains: header, APPROX\_WHITE\_POINT, COLOR_REP (iRGB), expected-values flag, fields, and patch data.

### 2. TI2 (.ti2)
Contains: header, APPROX\_WHITE\_POINT, strip/patch index patterns, STEPS\_IN\_PASS, PASSES\_IN\_STRIPS2, index order, and full color data including Lab.

`SAMPLE_LOC` uses labels like “A1”, “D14”, etc.

### 3. CSV (.csv)
Space-separated table with SAMPLE_ID, SAMPLE_LOC, and selected color values.

# Patch Value Ranges
RGB: 0–100
XYZ: typically 0–100
Lab_L: 0–100
Lab_a/b: approx. –128 to +128

All values have six decimal places.

# Valid Color Blocks
`--output_color_space` must include at least one of: RGB, XYZ, LAB

Example:
`--output_color_space RGB,XYZ,LAB`

Order controls column order in TI1/TI2/CSV.

# Dependencies
- Python 3.8+
- numpy
- Pillow
- colormath

# Command-line Arguments
| Argument | Description |
|---------|-------------|
| `--image` / `-i` | Path to input image containing colour patch grid. Input image must be a display-referred D65-based RGB image encoded with a 2.2 gamma transfer curve. |
| `--image_color_space` | **[Optional]** Input image device colour space; choices: `srgb` (default), `adobergb`. Determines the RGB→XYZ conversion matrices used by colormath. |
| `--patch_first_xy` | `"X,Y"` coordinates (floats or ints) of the **centre** of the first patch (top-left patch). Origin is top-left of the image. |
| `--patch_last_xy` | `"X,Y"` coordinates (floats or ints) of the **centre** of the last patch (bottom-right patch). |
| `--patch_width_height_ratio` | Width ÷ height (W/H) ratio of a single patch. Examples: `1.0` → square patches; `1.4` → width is 1.4× height. |
| `--num_cols` | Number of columns in the grid. |
| `--num_rows` | Number of rows in the grid. |
| `--row_labels`, `--col_labels` | Define the label sequences. Must match `num_rows` / `num_cols` exactly. Supported formats:<br> **Numeric:** `1-27`, `03-15`, `0001-0120` (zero-padding preserved).<br> **Alphabetic:** `A-Z`, `A-AA`, `BQ-CF`, up to `ZZ`.<br> Notes: Alphabetic and numeric are mutually exclusive between row and column—if rows use alphabetic labels, columns must be numeric, and vice-versa. |
| `--patch_label_order` | Determines how the patch label is formed: `col_then_row` → (col)(row), `row_then_col` → (row)(col). |
| `--output_color_space` | Comma-separated list. Allowed tokens: `RGB`, `XYZ`, `LAB`. Specifies which colour spaces to output in the generated file, and in what order. TI1 includes only RGB and/or XYZ; TI2/CSV includes any defined token.<br> Examples: `RGB,XYZ` → TI1 includes RGB,XYZ; TI2/CSV includes RGB,XYZ. `XYZ,LAB,RGB` → TI1 includes XYZ,RGB; TI2/CSV includes XYZ,LAB,RGB.<br> Output columns follow the specified sequence. |
| `--sample_fraction` | Fraction of patch to sample; default **0.20** (20%). Sampling square is centered on patch centre and clamped to at least **3×3 pixels** and at most **60%** of patch size. |
| `--sample_mode` | `"mean"` (default) or `"median"` for robust sampling. |
| `--output` | **[Optional]** Output filename base (defaults to `<imagebasename>.ti1/.ti2/.csv`). |

# Notes on Accuracy
- Uses 16-bit RGB internally.
- 8-bit images get scaled (value × 257).
- Color conversions use double precision.
- XYZ scaled so Yn=100.

# Patch Sampling Method
1. Compute patch center via interpolation.
2. Define sampling region size (min 3px, max 60% patch).
3. Extract pixels.
4. Average RGB in 16-bit domain.
5. Convert using colormath (RGB→XYZ D50 → Lab D50).

# Output File Format (TI2)
```
CTI2

DESCRIPTOR "Argyll Calibration Target chart information 2"
ORIGINATOR "read_image_patch_colors.py"
CREATED "<timestamp>"
ACCURATE_EXPECTED_VALUES "true"
APPROX_WHITE_POINT "<Xw Yw Zw>"
COLOR_REP "iRGB"
STEPS_IN_PASS "<ROWS>"
PASSES_IN_STRIPS2 "<COLS>"
STRIP_INDEX_PATTERN "<row_pattern>"
PATCH_INDEX_PATTERN "<column_pattern>"
INDEX_ORDER "STRIP_THEN_PATCH"

NUMBER_OF_FIELDS <AAA>
BEGIN_DATA_FORMAT
SAMPLE_ID SAMPLE_LOC [colour fields...]
END_DATA_FORMAT

NUMBER_OF_SETS <BBB>
BEGIN_DATA
  <patch_index> "<label>" <RGB/XYZ/Lab values...>
  ...
END_DATA
```

Fields appear in the order defined by `output_color_space`.

Examples:

- `RGB,XYZ` → RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z
- `LAB` → LAB_L LAB_A LAB_B

# Error Conditions
Script stops if:

- Label mismatch
- Invalid color tokens
- sample_fraction invalid
- Grid geometry cannot be computed
- File write error

# Workflow Summary
1. Validate dependencies
2. Parse CLI args
3. Load image
4. Compute grid geometry
5. Generate labels
6. Sample patches
7. Convert RGB→XYZ→Lab
8. Compute white/black points
9. Generate metadata
10. Write TI1/TI2/CSV

# Examples
## 9×12 SpyderPrint target with numeric rows, alphabetic columns
```
python3 read_image_patch_colors.py \
  --image "EZ 729 Colors Plus Grays 1 of 9 (108-patches).tif" \
  --patch_first_xy 111,64 \
  --patch_last_xy 1039,746 \
  --patch_width_height_ratio 1.89286 \
  --num_cols 9 \
  --num_rows 12 \
  --row_labels 1-12 \
  --col_labels A-I \
  --patch_label_order col_then_row \
  --output_color_space RGB,XYZ
```
## 10×12 SpyderPrint target with numeric rows, alphabetic columns
```
python3 read_image_patch_colors.py \
  --image "EZ 729 Colors Plus Grays 8 of 9 (120-patches).tif" \
  --patch_first_xy 87,67 \
  --patch_last_xy 1023,749 \
  --patch_width_height_ratio 1.7142857 \
  --num_cols 10 \
  --num_rows 12 \
  --row_labels 1-12 \
  --col_labels A-J \
  --patch_label_order col_then_row \
  --output_color_space RGB,XYZ
```
## 15×15 SpyderPrint target with numeric rows, alphabetic columns
```bash
python3 read_image_patch_colors.py \
  --image "High Quality Target (225-patches).tif" \
  --patch_first_xy 61,50 \
  --patch_last_xy 848,612 \
  --patch_width_height_ratio 1.37837838 \
  --num_cols 15 \
  --num_rows 15 \
  --row_labels 1-15 \
  --col_labels A-O \
  --patch_label_order col_then_row \
  --output_color_space RGB,XYZ
```
## 17×14 SpyderPrint target with numeric rows, alphabetic columns
```bash
python3 read_image_patch_colors.py \
  --image "Expert Target Plus Grays 4 of 4 (238-patches).tif" \
  --patch_first_xy 56,51 \
  --patch_last_xy 847,592 \
  --patch_width_height_ratio 1.09803921568 \
  --num_cols 17 \
  --num_rows 14 \
  --row_labels 1-14 \
  --col_labels A-Q \
  --patch_label_order col_then_row \
  --output_color_space RGB,XYZ
```

## 18×14 SpyderPrint target with numeric rows, alphabetic columns
**Note!**<br>
This script requires specification of first and last patch of a full grid. Thus, last 9 patches of last row must be manually removed from files as they do not exist in image.

```bash
python3 read_image_patch_colors.py \
  --image "Expert Target Page 1 of 3 (243-patches).tif" \
  --patch_first_xy 48,45 \
  --patch_last_xy 872,604 \
  --patch_width_height_ratio 1.075 \
  --num_cols 18 \
  --num_rows 14 \
  --row_labels 1-14 \
  --col_labels A-R \
  --patch_label_order col_then_row \
  --output_color_space RGB,XYZ
```
## 27×27 SpyderPrint target with numeric rows, alphabetic columns
**Note!**<br>
This image uses a special character on last column label, which must be manually changed in the ti2 file.
```
python3 read_image_patch_colors.py \
  --image "Expert Target (large)(729-patches).tif" \
  --patch_first_xy 56,38 \
  --patch_last_xy 1234,921 \
  --patch_width_height_ratio 1.4 \
  --num_cols 27 \
  --num_rows 27 \
  --row_labels 1-27 \
  --col_labels A-AA \
  --patch_label_order col_then_row \
  --output_color_space RGB,XYZ
```
## 36×26 LaserSoft target with numeric rows, alphabetic columns
```
python3 read_image_patch_colors.py \
  --image "LaserSoft_Advanced_Target-display.tif" \
  --patch_first_xy 168,168 \
  --patch_last_xy 1335,935 \
  --patch_width_height_ratio 1.0 \
  --num_cols 36 \
  --num_rows 24 \
  --row_labels A-X \
  --col_labels 01-36 \
  --patch_label_order row_then_col \
  --output_color_space RGB,XYZ
```
