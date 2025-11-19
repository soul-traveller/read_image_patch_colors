#!/usr/bin/env python3
"""
read_image_patch_colors.py
Version: 1.4.2

================================================================================
READ_IMAGE_PATCH_COLORS.PY — COMPLETE DOCUMENTATION
================================================================================

OVERVIEW
--------
This program extracts color values from a rectangular grid of color patches in
an image, computes colorimetric values (RGB percentages, XYZ, Lab), applies row/
column labeling rules, and writes three output files:

    • ArgyllCMS .ti1 file
    • ArgyllCMS .ti2 file
    • CSV file (space-separated)

It supports sampling in mean or median mode, configurable patch geometry,
numeric or alphabetic labeling patterns, and RGB/XYZ/Lab output combinations.

The script is designed for use in color-management workflows where printed color
targets must be scanned or photographed and converted into Argyll measurement
files.

--------------------------------------------------------------------------------
IMAGE INPUT
-----------
Accepted image types: any PIL-compatible raster file
Color depth: handled as 16-bit RGB internally
Color space:
    • "srgb"      (default)
    • "adobergb"

Parameter: --image / -i
Provides the path to the patch-grid image.

--------------------------------------------------------------------------------
GRID GEOMETRY
-------------
Geometry is computed from:
    • --patch_first_xy       top-left patch center (x,y)
    • --patch_last_xy        bottom-right patch center (x,y)
    • --num_cols             number of columns
    • --num_rows             number of rows
    • --patch_width_height_ratio  width/height ratio of a patch

The script calculates patch center coordinates and patch boundaries for sampling.

--------------------------------------------------------------------------------
SAMPLING
--------
Parameter: --sample_fraction
Range:     0 < f ≤ 0.6
Defines the fraction of the *smaller* patch dimension used for sampling.

Sampling modes:
    • "mean"
    • "median"

The sampling area is square and centered within each patch.

--------------------------------------------------------------------------------
LABELING RULES
--------------
Row and column labels are generated from:
    • --row_labels    e.g. "A-D" or "1-20"
    • --col_labels    e.g. "A-K" or "1-32"

Accepted patterns:
    • Alphabetic:  "A", "B", … "Z", "AA", "AB", …
    • Numeric:     integers ≥ 0

Patch label order is controlled by:
    --patch_label_order  = "col_then_row" or "row_then_col"

Row/column label counts must match grid dimensions exactly.

--------------------------------------------------------------------------------
COLOR PROCESSING
----------------
Each sampled RGB triplet is first normalized to percent:
    RGB_percent = (R_16bit / 65535) * 100

Valid range:  0.0 – 100.0

Color conversions (XYZ, Lab) are computed through colormath using the chosen
source color space ("srgb" or "adobergb").

Ranges:
    XYZ_X, XYZ_Y, XYZ_Z:  unbounded positive; typical printed values 0–100
    Lab_L:                0–100
    Lab_a, Lab_b:         approx. -128 to +128

The patch with the highest Y (XYZ_Y) is stored as APPROX_WHITE_POINT.

--------------------------------------------------------------------------------
OUTPUT FILES
------------
Three files are produced:

1. TI1 FILE (.ti1)
------------------
Contains:
    • Header (descriptor, originator, timestamp)
    • APPROX_WHITE_POINT
    • COLOR_REP (always "iRGB")
    • Accurate expected values flag
    • Data fields (SAMPLE_ID, RGB, optionally XYZ)
    • Patch data in Argyll .ti1 format

2. TI2 FILE (.ti2)
------------------
Contains:
    • Header (descriptor, originator, timestamp)
    • APPROX_WHITE_POINT
    • STRIP and PATCH index patterns inferred from row/col labels
    • STEPS_IN_PASS     = num_rows
    • PASSES_IN_STRIPS2 = num_cols
    • Index order       = STRIP_THEN_PATCH
    • Data fields (SAMPLE_ID, SAMPLE_LOC, RGB, XYZ, Lab)
    • Complete patch data

SAMPLE_LOC corresponds to labels such as "A1", "D14", "C07", etc.

3. CSV FILE (.csv)
------------------
A space-separated ASCII data table containing:
    SAMPLE_ID, SAMPLE_LOC, and selected color values.

--------------------------------------------------------------------------------
PATCH VALUE RANGES
------------------
RGB_R/G/B:   0.000000 – 100.000000
XYZ_X/Y/Z:   Typically 0–100 but depends on color space
LAB_L:       0.000000 – 100.000000
LAB_A/B:     Approximately -128 – +128

Values are always written with 6 decimal places.

--------------------------------------------------------------------------------
VALID COLOR BLOCKS
------------------
--output_color_space must be a comma-separated list including at least one of:
    RGB, XYZ, LAB

Example:
    --output_color_space RGB,XYZ,LAB

The same set controls which fields, and in what order, that appear in TI1, TI2, and CSV.

--------------------------------------------------------------------------------
DEPENDENCIES
------------
The script requires:

    • Python 3.8+
    • numpy
    • Pillow (PIL)
    • colormath

If missing, the script prints installation instructions and exits.

--------------------------------------------------------------------------------
ERROR CONDITIONS
----------------
The script stops with a readable message when:
    • Labels do not match row/column counts
    • Output color space tokens are invalid
    • sample_fraction exceeds allowed range
    • Grid geometry cannot be computed
    • Any I/O error occurs during file writing

--------------------------------------------------------------------------------
WORKFLOW SUMMARY
----------------
1. Validate dependencies
2. Parse CLI arguments
3. Load image as 16-bit RGB
4. Compute patch grid geometry
5. Generate label sets
6. Sample each patch
7. Convert sampled RGB to XYZ and Lab
8. Compute white/black points
9. Generate patch metadata
10. Write TI1, TI2, and CSV output files

================================================================================
COMMAND-LINE ARGUMENTS
================================================================================
--image / -i                    Path to input image containing colour patch grid.
                                Input image must be a display ready D65 image with gamma 2.2 applied.
--image_color_space             [Optional] input image device colour space; choices:
                                srgb (default), adobergb.
                                Determines the RGB->XYZ conversion matrices used by colormath.
--patch_first_xy                "X,Y" coordinates (floats or ints) of the ***centre***
                                of the first patch (top-left patch). Origin is top-left of the image.
--patch_last_xy                 "X,Y" coordinates (floats or ints) of the ***centre*** of the last patch (bottom-right patch).
--patch_width_height_ratio      Width ÷ height (W/H) ratio of a single patch.
                                    Examples:
                                      1.0   → square patches
                                      1.4   → width is 1.4× height
--num_cols                      Number of columns in the grid.
--num_rows                      Number of rows in the grid.
--row_labels
--col_labels
                                Define the label sequences. Must match num_rows / num_cols exactly.
                                Supported formats:
                                • Numeric: 1-27, 03-15, 0001-0120 (zero-padding preserved)
                                • Alphabetic: A-Z, A-AA, BQ-CF, up to ZZ
                                Notes:
                                • Alphabetic and numeric are mutually exclusive between row and column:
                                  If rows use alphabetic labels, columns must be numeric, and vice-versa.
--patch_label_order             Determines how the patch label is formed:
                                  col_then_row → "<col><row>"
                                  row_then_col → "<row><col>"
--output_color_space           Comma-separated list. Allowed tokens: RGB, XYZ, LAB
                                Specifies **which colour spaces** to output in generated file, and in what order.
                                TI1 only includes RGB and/or XYZ. TI2/CSV includes any defined token.
                                Example tokens:
                                  RGB,XYZ       →   TI1 includes RGB,XYZ; TI2/CSV includes RGB,XYZ
                                  XYZ,LAB,RGB   →   TI1 includes XYZ,RGB; TI2/CSV includes XYZ,LAB,RGB
                                Output columns follow the specified sequence.
--sample_fraction               Fraction of patch to sample; default 0.20 (20%).
                                The sampling square is centred on the patch centre and is clamped to
                                at least 3×3 pixels and at most 60% of patch size.
--sample_mode                    "mean" (default) or "median" for robust sampling
--output                         [Optional] output filename (defaults to <imagebasename>.ti1/ti2/csv)

================================================================================
COLOR / SCALING CONVENTIONS IN OUTPUT FILES (Argyll-compatible)
================================================================================
- RGB are written as **iRGB** percentage values: 0.000000 .. 100.000000
  (Argyll convention; device nominal percentages).

- CIEXYZ in output:
  * CIEXYZ (CIE 1931 2° standard observer)
  * **D50** white reference (Profile Connection Space in ICC)
  * **Scaled so that white has Y = 100.000000**
  * Numeric precision: six decimal places.

- Lab in output:
  * CIE L*a*b* (1976), computed using **D50** reference white (Xn, Yn, Zn)
  * L in range ~0..100, a and b around typical -128..128 ranges.
  * Numeric precision: six decimal places.

- All colorimetric output is **D50** and **linear** (no gamma applied in XYZ or Lab).
  The conversion pipeline in this script applies TRC/matrix conversions according to
  the selected `--image_color_space` using python-colormath. Chromatic adaptation to
  D50 is performed by colormath during conversion when required.

================================================================================
NOTES ON ACCURACY
================================================================================
- This script seeks to minimize rounding/quantization error:
  * If the image has 16-bit/channel data it will be used directly.
  * If the image is 8-bit/channel (most typical images), it will be scaled up
    exactly to 16-bit (value * 257) before colour conversion to reduce rounding
    / mapping error in conversions.
  * Conversion from RGB -> XYZ -> Lab is done using python-colormath with double
    precision floats; XYZ is scaled so Yn = 100 and written with 6 decimal places.

- This script produces **colorimetrically correct** PCS values
  for the image that was measured.

================================================================================
Patch Sampling Method
================================================================================
For each patch:

1. Patch centre (Cx, Cy) is computed by linear interpolation between
   patch_first_xy → patch_last_xy across the grid.
2. A square sampling region of size `sample_size × sample_size`
   (odd number, min 3, max 60% of patch) is centred at (Cx, Cy).
3. The script extracts all pixels inside this region.
4. The sampled RGB values are averaged in 16-bit integer space (8-bit inputs are scaled to 16-bit via *257).
5. Colour conversions are applied using python-colormath:
      averaged 16-bit RGB (scaled to 0..1 via /65535)
          -> colormath: RGB -> CIEXYZ (D50), then XYZ -> Lab (D50)

================================================================================
Output File Format (ArgyllCMS TI2)
================================================================================
Output files:
    <input_filename>.ti1
    <input_filename>.ti2
    <input_filename>.csv

Structure:

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

Where:
- SAMPLE_ID is the numeric row-major patch index starting at 1
- SAMPLE_LOC is the label, e.g. "A5" or "15C"
- AAA = number of fields (depends on output_color_space)
- BBB = number of patches (num_rows × num_cols)
- ROWS = number of rows (num_rows)
- COLS = number of columns (num_cols)

Colour fields appear **in the same order specified by output_color_space**.
Examples:

    output_color_space="RGB,XYZ"
    → RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z

    output_color_space="LAB"
    → LAB_L LAB_A LAB_B

================================================================================
Examples
================================================================================

Image with 9 by 12 grid chart, alphabetic columns and numeric rows:
    python3 read_image_patch_colors.py \
      --image "EZ 729 Colors Plus Grays 1 of 9.tif" \
      --patch_first_xy 111,64 \
      --patch_last_xy 1039,746 \
      --patch_width_height_ratio 1.89286 \
      --num_cols 9 \
      --num_rows 12
      --row_labels 1-12 \
      --col_labels A-I \
      --patch_label_order col_then_row \
      --output_color_space RGB,XYZ

This will write:
 - "EZ 729 Colors Plus Grays 1 of 9.ti2"
 in the current directory.

Image with 10 by 12 grid chart, alphabetic columns and numeric rows:
     python3 read_image_patch_colors.py \
      --image "EZ 729 Colors Plus Grays 8 of 9.tif" \
      --patch_first_xy 87,67 \
      --patch_last_xy 1023,749 \
      --patch_width_height_ratio 1.7142857 \
      --num_cols 10 \
      --num_rows 12
      --row_labels 1-12 \
      --col_labels A-J \
      --patch_label_order col_then_row \
      --output_color_space RGB,XYZ

Image with 27 by 27 grid chart, alphabetic columns and numeric rows:
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
END OF HEADER
"""

from __future__ import annotations
import argparse
import math
import time
import os
import sys
import re
import platform
from dataclasses import dataclass
from typing import List, Tuple, Sequence, Optional

# ---------------- Preflight ----------------
# Check that required Python packages are installed.
# Input: None
# Output: list of missing package names (list[str]).
# Packages checked: colormath, numpy, Pillow
def check_full_packages() -> list[str]:
    missing = []
    import importlib.util

    if importlib.util.find_spec("colormath") is None:
        missing.append("colormath")
    if importlib.util.find_spec("numpy") is None:
        missing.append("numpy")
    if importlib.util.find_spec("PIL") is None:
        missing.append("Pillow")
    return missing


# Print platform-specific pip install instructions for missing packages.
# Input: packages (list[str]) — names of missing packages
# Output: prints instructions to stdout; no return value.
def print_install_instructions(packages: list[str]) -> None:
    system = platform.system()
    print(f"Missing packages detected: {', '.join(packages)}")
    print("Install them manually with pip, for example:")
    if system == "Windows":
        print(f"  python -m pip install {' '.join(packages)}")
    else:
        print(f"  python3 -m pip install {' '.join(packages)}")


missing_packages = check_full_packages()
if missing_packages:
    print_install_instructions(missing_packages)
    sys.exit(1)


# ---------------- Imports AFTER preflight ----------------
from PIL import Image
import numpy as np
from colormath.color_objects import sRGBColor, AdobeRGBColor, XYZColor, LabColor
from colormath.color_conversions import convert_color

# ---------- Constants ----------
DEBUG = False
DEBUG_PRINT_LIMIT = 10
DEBUG_AVG_COUNTER = 0
DEBUG_SAMPLE_COUNTER = 0
# Reference white D50 (used in Lab conversion)
D50_X = 96.422
D50_Y = 100.000
D50_Z = 82.521
DEFAULT_SAMPLE_FRACTION = 0.20
EPS = (6.0 / 29.0) ** 3   # Threshold for linearization in Lab
K = 24389.0 / 27.0        # Linear coefficient in Lab conversion


# ---------- Utilities ----------
# Create a formatted timestamp string suitable for TI files (CTI1/CTI2).
# Input: t_struct (time.struct_time) — typically time.localtime().
# Output: string like 'Wed Mar 20 14:30:05 2024'.
def make_created_timestamp(t_struct: time.struct_time) -> str:
    wday = time.strftime("%a", t_struct)
    mon = time.strftime("%b", t_struct)
    day = int(time.strftime("%d", t_struct))
    hhmmss = time.strftime("%H:%M:%S", t_struct)
    year = time.strftime("%Y", t_struct)
    day_str = f"{day:2d}"
    return f'{wday} {mon} {day_str} {hhmmss} {year}'


# Parse a string containing "X,Y" into floats.
# Input: s (str) e.g. "123.4,567.8".
# Output: tuple (float, float) representing (X, Y).
# Raises argparse.ArgumentTypeError on invalid input.
def parse_xy(s: str) -> Tuple[float, float]:
    try:
        parts = s.split(',')
        if len(parts) != 2:
            raise ValueError()
        return float(parts[0].strip()), float(parts[1].strip())
    except Exception:
        raise argparse.ArgumentTypeError("XY must be X,Y")


# Clamp an integer value to the inclusive range [amin, amax].
# Input: v (int), amin (int), amax (int)
# Output: int within range.
def clamp_int(v: int, amin: int, amax: int) -> int:
    return max(amin, min(amax, v))


# ---------- Label handling ----------
# Convert alphabetic Excel-like label (A..Z, AA..ZZ) to 1-based index.
# Input: s (str) uppercase letters A..ZZ
# Output: int index where A=1, B=2, ..., ZZ=702
def alpha_to_int(s: str) -> int:
    v = 0
    for c in s:
        v = v * 26 + (ord(c) - ord('A') + 1)
    return v


# Convert 1-based integer index to alphabetic label (A..ZZ).
# Input: i (int) in 1..702
# Output: string label ("A", "Z", "AA", ...)
# Raises ValueError if out of range.
def int_to_alpha(i: int) -> str:
    if not (1 <= i <= 26*27):
        raise ValueError("Alphabetic label index out of range (1–702)")
    out = ""
    while i > 0:
        i -= 1
        out = chr((i % 26) + ord('A')) + out
        i //= 26
    return out


# Parse a label range specification into a list of labels and pad length.
# Input: spec (str) e.g. "1-20", "03-15", "A-D".
# Output: (labels: list[str], pad: Optional[int]) where pad is int for numeric padding or None.
# Raises ValueError on invalid formats or reversed ranges.
def parse_label_range(spec: str) -> tuple[list[str], Optional[int]]:
    if "-" not in spec:
        raise ValueError("Label range must contain '-'")
    start, end = spec.split("-", 1)
    start, end = start.strip(), end.strip()

    # Numeric range (preserve zero padding)
    if re.fullmatch(r"\d+", start) and re.fullmatch(r"\d+", end):
        pad = len(start)
        s, e = int(start), int(end)
        if e < s:
            raise ValueError("Numeric label range end < start")
        return [str(i).zfill(pad) for i in range(s, e+1)], pad

    # Alphabetic range
    if re.fullmatch(r"[A-Z]+", start) and re.fullmatch(r"[A-Z]+", end):
        s, e = alpha_to_int(start), alpha_to_int(end)
        if e < s:
            raise ValueError("Alphabetic label range end < start")
        if e > 26*27:
            raise ValueError("Max alphabetic label is ZZ (702)")
        return [int_to_alpha(i) for i in range(s, e+1)], None

    raise ValueError(f"Invalid label range format: {spec}")


# ---------- Lab <-> XYZ helpers ----------
# Forward function f(t) used in Lab conversion.
# Input: t_arr (np.ndarray) ratios X/Xn, Y/Yn, Z/Zn (float)
# Output: f(t) as np.ndarray same shape — cubic root or linear approx.
def f_for_lab_arr(t_arr: np.ndarray) -> np.ndarray:
    t_arr = np.asarray(t_arr, dtype=np.float64)
    mask = t_arr > EPS
    out = np.empty_like(t_arr, dtype=np.float64)
    # Non-linear (cubic root) branch for t > EPS
    out[mask] = np.cbrt(t_arr[mask])
    # Linear approximation branch for small t
    out[~mask] = (K * t_arr[~mask] + 16.0) / 116.0
    return out


# Convert CIEXYZ to CIE Lab (D50 reference).
# Input: xyz (np.ndarray, shape (...,3)), values typically 0..100 for X,Y,Z (Y scaled to 100).
# Output: Lab array same shape: L in ~[0..100], a,b roughly [-128..128].
# Assumes reference white D50 (D50_X, D50_Y, D50_Z).
def xyz_to_lab_d50(xyz: np.ndarray) -> np.ndarray:
    xyz = np.asarray(xyz, dtype=np.float64)

    # Normalize by the D50 whitepoint (ratios)
    fx = f_for_lab_arr(xyz[..., 0] / D50_X)
    fy = f_for_lab_arr(xyz[..., 1] / D50_Y)
    fz = f_for_lab_arr(xyz[..., 2] / D50_Z)

    # Compute L*, a*, b*
    L = (116.0 * fy) - 16.0
    a = 500.0 * (fx - fy)
    b = 200.0 * (fy - fz)

    return np.stack([L, a, b], axis=-1)


# Inverse function used for Lab -> XYZ conversion.
# Input: f_arr (np.ndarray) f-values
# Output: t-values (np.ndarray) suitable for XYZ reconstruction.
def finv_arr(f_arr: np.ndarray) -> np.ndarray:
    f_arr = np.asarray(f_arr, dtype=np.float64)
    mask = f_arr > (6.0 / 29.0)
    out = np.empty_like(f_arr, dtype=np.float64)

    # Cubic inverse for larger f
    out[mask] = f_arr[mask] ** 3
    # Linear inverse for small f
    out[~mask] = (116.0 * f_arr[~mask] - 16.0) / K
    return out


# Convert Lab (D50) to CIEXYZ (scaled so Y = 100).
# Input: lab (np.ndarray, shape (...,3)) with L in 0..100
# Output: xyz (np.ndarray, same shape) with X,Y,Z scaled so Y=100.
def lab_to_xyz_d50(lab: np.ndarray) -> np.ndarray:
    lab = np.asarray(lab, dtype=np.float64)
    L = lab[..., 0]
    a = lab[..., 1]
    b = lab[..., 2]

    fy = (L + 16.0) / 116.0
    fx = fy + (a / 500.0)
    fz = fy - (b / 200.0)

    xr = finv_arr(np.asarray(fx))
    yr = finv_arr(np.asarray(fy))
    zr = finv_arr(np.asarray(fz))

    x = xr * (D50_X / 100.0)
    y = yr * (D50_Y / 100.0)
    z = zr * (D50_Z / 100.0)

    xyz100 = np.stack([x * 100.0, y * 100.0, z * 100.0], axis=-1)
    return xyz100


# ---------- Image I/O and averaging ----------
# Load an image and return a 16-bit RGB numpy array (HxWx3, uint16).
# Input: img_path (str)
# Output: tuple (arr16: np.ndarray HxWx3 uint16, bit_depth: int 8|16)
def load_image_as_16bit_rgb(img_path: str) -> Tuple[np.ndarray, int]:
    try:
        im = Image.open(img_path)
    except Exception as e:
        print('Error: Failed to open image:', e)
        sys.exit(1)

    # Convert to RGB (drops alpha, converts grayscale to RGB etc.)
    im = im.convert('RGB')
    arr = np.array(im)

    if arr.dtype == np.uint8:
        # Expand 8-bit to full 16-bit range exactly (value * 257)
        arr16 = (arr.astype(np.uint32) * 257).astype(np.uint16)
        return arr16, 8

    elif arr.dtype == np.uint16:
        # Already 16-bit per channel
        return arr.astype(np.uint16), 16

    else:
        # Floating or other range: clamp to 0..1 then scale
        arrf = np.clip(arr.astype(np.float64), 0.0, 1.0)
        arr16 = (arrf * 65535.0).round().astype(np.uint16)
        return arr16, 16


# Compute the RGB value of a patch using different aggregation modes.
# Input:
#    block : (H x W x 3) ndarray, 16-bit values already loaded as uint16
#    mode  : 'mean'   – simple arithmetic mean (not robust to outliers)
#            'median' – statistically robust, but less efficient for smooth gradients
#            'mad'    – MAD-based sigma-clipped robust mean (default)
# Output:
#    1x3 ndarray of float64 RGB values in the 0..65535 range
#
# Notes:
#   - 'mad' mode is recommended for scanned/photographic data where dust,
#     specks, and sensor defects introduce heavy-tailed noise.
#   - Only the selected mode is returned. There is no post-override.
def average_block_rgb16(block: np.ndarray, mode: str = 'mad') -> np.ndarray:
    global DEBUG_AVG_COUNTER

    if block.size == 0:
        return np.array([0.0, 0.0, 0.0], dtype=np.float64)

    arr = block.reshape(-1, 3).astype(np.float64)

    # ----- simple mean -----
    if mode == 'mean':
        return arr.mean(axis=0)

    # ----- simple median -----
    if mode == 'median':
        return np.median(arr, axis=0)

    # ----- robust MAD-sigma-clipped mode -----
    def robust_channel_mean(values: np.ndarray, k: float = 3.0, min_sigma: float = 1.0) -> float:
        if values.size == 0:
            return 0.0
        med = np.median(values)
        mad = np.median(np.abs(values - med))
        # Convert MAD → robust sigma estimate with minimum floor
        sigma = max(1.4826 * mad, min_sigma)
        # sigma-based threshold
        abs_thresh = k * sigma
        mask = np.abs(values - med) <= abs_thresh
        if mask.sum() == 0:
            # If all rejected, use median
            return float(med)
        return float(values[mask].mean())

    # Apply per-channel MAD clipping
    r = robust_channel_mean(arr[:, 0])
    g = robust_channel_mean(arr[:, 1])
    b = robust_channel_mean(arr[:, 2])
    robust_mean = np.array([r, g, b], dtype=np.float64)

    # Debug native and robust mean for first 10 MAD patches
    if DEBUG and mode == 'mad' and DEBUG_AVG_COUNTER < DEBUG_PRINT_LIMIT:
        patch_num = DEBUG_AVG_COUNTER + 1
        print(f"--- PATCH {patch_num} ---")
        print(f"Naive mean (16-bit): {arr.mean(axis=0)}")
        print(f"Robust mean (16-bit): {robust_mean}")
        # Show MAD / sigma info per channel
        for i, ch in enumerate(["R", "G", "B"]):
            med = np.median(arr[:, i])
            mad = np.median(np.abs(arr[:, i] - med))
            sigma = max(1.4826 * mad, 1.0)
            mask = np.abs(arr[:, i] - med) <= 3 * sigma
            rejected = (~mask).sum()
            print(f"{ch}-channel: median={med:.1f}, MAD={mad:.1f}, sigma={sigma:.1f}, outliers removed={rejected}")
        DEBUG_AVG_COUNTER += 1

    return robust_mean


# Convert 16-bit RGB to Argyll iRGB percentage (0.0 .. 100.0)
# Input: rgb16 (3 values 0..65535)
# Output: tuple of floats in 0..100
def rgb16_to_argyll_percent(rgb16: Sequence[float]) -> Tuple[float, float, float]:
    factor = 100.0 / 65535.0
    return float(rgb16[0] * factor), float(rgb16[1] * factor), float(rgb16[2] * factor)


# Map input color space name to a colormath RGB color class and kwargs.
# Input: space_name (str) - 'srgb' or 'adobergb'
# Output: (ColorClass, kwargs dict)
def get_rgb_class_and_kwargs(space_name: str):
    s = space_name.lower()
    if s == 'srgb':
        return sRGBColor, {}
    if s in ('adobergb', 'adobergb1998', 'adobergb-1998'):
        return AdobeRGBColor, {}
    raise ValueError(f"Unsupported image_color_space: {space_name}")


# Convert averaged 16-bit RGB to CIEXYZ (D50, scaled so Y=100) and Lab (D50).
# Input: avg_rgb16 sequence of 3 values (0..65535), image_color_space (str)
# Output: (xyz100: ndarray(3,), lab_float: ndarray(3,))
def rgb16_to_xyz_lab(avg_rgb16: Sequence[float], image_color_space: str) -> Tuple[np.ndarray, np.ndarray]:
    # Scale to 0..1 for colormath
    rgb01 = np.asarray(avg_rgb16, dtype=np.float64) / 65535.0

    ColorClass, kwargs = get_rgb_class_and_kwargs(image_color_space)
    rgb_obj = ColorClass(rgb01[0], rgb01[1], rgb01[2], **kwargs)

    # Convert RGB -> XYZ, request D50 illuminant if supported
    try:
        xyz_obj = convert_color(rgb_obj, XYZColor, target_illuminant='d50')
    except Exception:
        xyz_obj = convert_color(rgb_obj, XYZColor)

    # Extract XYZ (colormath returns 0..1-ish for PCS; we scale to Y=100 below)
    X = float(getattr(xyz_obj, 'xyz_x', xyz_obj.xyz_x))
    Y = float(getattr(xyz_obj, 'xyz_y', xyz_obj.xyz_y))
    Z = float(getattr(xyz_obj, 'xyz_z', xyz_obj.xyz_z))

    # Scale to XYZ where Yn = 100
    xyz100 = np.array([X * 100.0, Y * 100.0, Z * 100.0], dtype=np.float64)

    # Convert XYZ -> Lab (D50)
    try:
        lab_obj = convert_color(xyz_obj, LabColor, target_illuminant='d50')
    except Exception:
        lab_obj = convert_color(xyz_obj, LabColor)

    L = float(getattr(lab_obj, 'lab_l', lab_obj.lab_l))
    a = float(getattr(lab_obj, 'lab_a', lab_obj.lab_a))
    b = float(getattr(lab_obj, 'lab_b', lab_obj.lab_b))

    lab_float = np.array([L, a, b], dtype=np.float64)

    return xyz100, lab_float


# ---------- Main Processing Data Classes ----------
@dataclass
class PatchInfo:
    """Container for a single patch measurement.

    Attributes:
        index: 1-based patch index in iteration order.
        label: human-readable label (e.g. "A1").
        rgb16: averaged 16-bit RGB values (numpy array, 0..65535).
        rgb_percent: tuple of R,G,B in Argyll percent (0..100).
        xyz100: CIEXYZ scaled so Y=100 (numpy array).
        lab: CIE Lab (numpy array).
    """
    index: int
    label: str
    rgb16: np.ndarray
    rgb_percent: Tuple[float, float, float]
    xyz100: np.ndarray
    lab: np.ndarray


@dataclass
class GridGeometry:
    """Geometry and sampling parameters for the patch grid.

    Attributes mirror the CLI inputs and derived sample size/spacing.
    """
    fx: float
    fy: float
    lx: float
    ly: float
    num_cols: int
    num_rows: int
    ratio: float
    sample_fraction: float
    patch_width: float
    patch_height: float
    sample_size: int
    half: int
    dx: float
    dy: float


# -------- Patch Grid Utilities --------
# Compute grid geometry and derived sampling size from CLI args.
# Input: args with attributes patch_first_xy, patch_last_xy, patch_width_height_ratio, num_cols, num_rows, sample_fraction.
# Output: GridGeometry instance with sample_size, half, dx, dy, patch sizes.
def compute_grid_geometry(args) -> GridGeometry:
    fx, fy = args.patch_first_xy
    lx, ly = args.patch_last_xy
    ratio = float(args.patch_width_height_ratio)
    num_cols, num_rows = int(args.num_cols), int(args.num_rows)
    sample_fraction = float(args.sample_fraction)

    # Validate basic ranges
    if num_cols < 2 or num_rows < 2:
        raise ValueError('num_cols and num_rows must both be >=2')
    if sample_fraction <= 0.0 or sample_fraction > 0.6:
        raise ValueError('sample_fraction must be in (0,0.6]')
    if ratio <= 0.0:
        raise ValueError('patch_width_height_ratio must be >0')

    # Compute center-to-center horizontal spacing
    patch_width = (lx - fx) / float(num_cols - 1)      # distance between patch centers
    patch_height = patch_width / ratio

    # Determine sampling square size (clamped between 3 and 60% of patch size)
    small_patch_side = min(patch_width, patch_height)
    sample_size = max(3, int(round(small_patch_side * sample_fraction)))
    max_allowed = int(math.floor(small_patch_side * 0.6))
    if max_allowed < 3:
        max_allowed = 3
    sample_size = min(sample_size, max_allowed)
    if sample_size % 2 == 0:
        sample_size += 1

    half = sample_size // 2
    dx = (lx - fx) / float(num_cols - 1)
    dy = (ly - fy) / float(num_rows - 1)

    # Friendly logging of geometry
    print(f'Computed patch_width = {patch_width:.6f} px')
    print(f'Computed patch_height = {patch_height:.6f} px')
    print(f'Sampling a central square of {sample_size}x{sample_size} pixels for each patch')

    return GridGeometry(
        fx=fx, fy=fy, lx=lx, ly=ly,
        num_cols=num_cols, num_rows=num_rows,
        ratio=ratio, sample_fraction=sample_fraction,
        patch_width=patch_width, patch_height=patch_height,
        sample_size=sample_size, half=half,
        dx=dx, dy=dy
    )


# Sample a centered square region from the image and aggregate pixels.
# Input: img16 (HxWx3 uint16 array), cx, cy center coords (float), half (int), sample_mode ('mean'|'median')
# Output: 1x3 ndarray averaged rgb16 values (float)
def sample_patch(img16: np.ndarray, cx: float, cy: float, half: int, sample_mode: str) -> np.ndarray:
    global DEBUG_SAMPLE_COUNTER

    height, width, _ = img16.shape
    # Compute pixel bounds of the sampling square
    x0 = int(round(cx)) - half
    x1 = int(round(cx)) + half + 1
    y0 = int(round(cy)) - half
    y1 = int(round(cy)) + half + 1

    # Clamp to image bounds
    x0c = clamp_int(x0, 0, width - 1)
    x1c = clamp_int(x1, 1, width)
    y0c = clamp_int(y0, 0, height - 1)
    y1c = clamp_int(y1, 1, height)

    if x0c >= x1c or y0c >= y1c:
        # Degenerate/empty block, return zeros
        return np.zeros(3, dtype=np.float64)

    block = img16[y0c:y1c, x0c:x1c, :]

    # ---- DEBUG BLOCK: prints for first 10 patches ----
    if DEBUG and DEBUG_SAMPLE_COUNTER < DEBUG_PRINT_LIMIT:
        patch_num = DEBUG_SAMPLE_COUNTER + 1
        flat = block.reshape(-1, 3).astype(np.float64)
        mins = flat.min(axis=0)
        maxs = flat.max(axis=0)
        medians = np.median(flat, axis=0)
        means = flat.mean(axis=0)
        stds = flat.std(axis=0)

        print(f"\n=== PATCH {patch_num} ===")
        print(f"Block shape: {block.shape}, dtype: {block.dtype}")
        print(f"Min / Max per channel: {mins} / {maxs}")
        print(f"Median / Mean / Std per channel: {medians} / {means} / {stds}")

        # Show a small random sample of pixel values
        sample_indices = np.random.choice(flat.shape[0], min(30, flat.shape[0]), replace=False)
        print("Random Pixel samples (R,G,B):", flat[sample_indices])

        try:
            from PIL import Image
            block8 = np.round(block / 257.0).clip(0, 255).astype(np.uint8)
            debug_path = f"debug_patch_{patch_num}.png"
            Image.fromarray(block8, mode='RGB').save(debug_path)
            print(f"Saved block image: {debug_path}")
        except Exception as e:
            print(f"Failed to save block image: {e}")

        DEBUG_SAMPLE_COUNTER += 1
    # ---- END DEBUG BLOCK ----

    avg_rgb16 = average_block_rgb16(block, mode=sample_mode)
    return avg_rgb16


# Extract measurements for all patches in the grid.
# Traverses the patch grid either row-major (default: row by row) or column-major (column by column),
# optionally reversing the order of patches per row (row-major) or per column (column-major) if mirror_output=True.
#
# Input:
#   img16             : HxWx3 uint16 array representing the image
#   geometry          : GridGeometry object with patch layout and spacing
#   row_labels        : list of row labels (strings)
#   col_labels        : list of column labels (strings)
#   row_pad           : optional int, whether row labels should be zero-padded
#   col_pad           : optional int, whether column labels should be zero-padded
#   patch_label_order : 'col_then_row' or 'row_then_col', determines label composition
#   image_color_space : string, color space of input image ('srgb', 'adobergb', etc.)
#   sample_mode       : 'mean' or 'median', aggregation method for patch sampling
#   output_order      : 'row_major' or 'column_major', controls traversal order
#   mirror_output     : bool, if True reverses order of patches per row (row-major) or per column (column-major)
#   rotate_grid       : integer, degrees to rotate grid clockwise.
#
# Output:
#   patches           : List[PatchInfo] objects with measured RGB, XYZ, LAB, label, and index
#   white_point_xyz   : XYZ coordinates of the patch with highest Y luminance
#   black_point_xyz   : XYZ coordinates of the patch with lowest Y luminance
def extract_patch_data(
    img16: np.ndarray,
    geometry: GridGeometry,
    row_labels: List[str],
    col_labels: List[str],
    row_pad: Optional[int],
    col_pad: Optional[int],
    patch_label_order: str,
    image_color_space: str,
    sample_mode: str,
    output_order: str,
    mirror_output: bool,
    rotate_grid: int
) -> Tuple[List[PatchInfo], np.ndarray, np.ndarray]:

    patches: List[PatchInfo] = []
    patch_idx = 0
    white_point_xyz = None
    black_point_xyz = None

    # ------------------------------------------------------
    # ROW-MAJOR ORDER
    # ------------------------------------------------------
    if output_order == "row_major":

        for row in range(geometry.num_rows):
            for col in range(geometry.num_cols):

                patch_idx += 1

                rlab = row_labels[row]
                clab = col_labels[col]
                if row_pad == 1:
                    rlab = str(int(rlab))
                if col_pad == 1:
                    clab = str(int(clab))
                patch_label = (
                    f"{clab}{rlab}"
                    if patch_label_order == "col_then_row"
                    else f"{rlab}{clab}"
                )

                cx = geometry.fx + col * geometry.dx
                cy = geometry.fy + row * geometry.dy

                avg_rgb16 = sample_patch(img16, cx, cy, geometry.half, sample_mode)
                rgb_percent = rgb16_to_argyll_percent(avg_rgb16)
                xyz100, lab_float = rgb16_to_xyz_lab(avg_rgb16, image_color_space)

                Y_lum = float(xyz100[1])
                if white_point_xyz is None or Y_lum > float(white_point_xyz[1]):
                    white_point_xyz = xyz100.copy()
                if black_point_xyz is None or Y_lum < float(black_point_xyz[1]):
                    black_point_xyz = xyz100.copy()

                patches.append(
                    PatchInfo(
                        index=patch_idx,
                        label=patch_label,
                        rgb16=avg_rgb16,
                        rgb_percent=rgb_percent,
                        xyz100=xyz100,
                        lab=lab_float
                    )
                )

        # ----- Rotate grid before mirroring -----
        if rotate_grid != 0:
            # After building patches, reshape into 2D grid
            patch_grid = np.array(patches, dtype=object).reshape((geometry.num_rows, geometry.num_cols))
            # Apply rotation
            k = rotate_grid // 90  # np.rot90 rotates counter-clockwise
            patch_grid = np.rot90(patch_grid, k=k)
            # Flatten back to list in row-major order
            patches = patch_grid.flatten().tolist()
            # Update indices
            for i, p in enumerate(patches, start=1):
                p.index = i

        # ------------------------------------------------------
        # Row-based mirroring
        # ------------------------------------------------------
        if mirror_output:
            mirrored = []
            for r in range(geometry.num_rows):
                start = r * geometry.num_cols
                end = start + geometry.num_cols
                row_patches = patches[start:end]
                mirrored.extend(row_patches[::-1])
            patches = mirrored

            for i, p in enumerate(patches, start=1):
                p.index = i

    # ------------------------------------------------------
    # COLUMN-MAJOR ORDER
    # ------------------------------------------------------
    else:  # column_major

        for col in range(geometry.num_cols):
            for row in range(geometry.num_rows):

                patch_idx += 1

                rlab = row_labels[row]
                clab = col_labels[col]
                if row_pad == 1:
                    rlab = str(int(rlab))
                if col_pad == 1:
                    clab = str(int(clab))
                patch_label = (
                    f"{clab}{rlab}"
                    if patch_label_order == "col_then_row"
                    else f"{rlab}{clab}"
                )

                cx = geometry.fx + col * geometry.dx
                cy = geometry.fy + row * geometry.dy

                avg_rgb16 = sample_patch(img16, cx, cy, geometry.half, sample_mode)
                rgb_percent = rgb16_to_argyll_percent(avg_rgb16)
                xyz100, lab_float = rgb16_to_xyz_lab(avg_rgb16, image_color_space)

                Y_lum = float(xyz100[1])
                if white_point_xyz is None or Y_lum > float(white_point_xyz[1]):
                    white_point_xyz = xyz100.copy()
                if black_point_xyz is None or Y_lum < float(black_point_xyz[1]):
                    black_point_xyz = xyz100.copy()

                patches.append(
                    PatchInfo(
                        index=patch_idx,
                        label=patch_label,
                        rgb16=avg_rgb16,
                        rgb_percent=rgb_percent,
                        xyz100=xyz100,
                        lab=lab_float
                    )
                )

        # ----- Rotate grid before mirroring -----
        if rotate_grid != 0:
            # After building patches, reshape into 2D grid
            patch_grid = np.array(patches, dtype=object).reshape((geometry.num_rows, geometry.num_cols))
            # Apply rotation
            k = rotate_grid // 90  # np.rot90 rotates counter-clockwise
            patch_grid = np.rot90(patch_grid, k=k)
            # Flatten back to list in column-major order
            patches = patch_grid.flatten().tolist()
            # Update indices
            for i, p in enumerate(patches, start=1):
                p.index = i

        # ------------------------------------------------------
        # Column-major mirroring (per column)
        # ------------------------------------------------------
        if mirror_output:
            mirrored = []
            for c in range(geometry.num_cols):
                col_indices = range(c * geometry.num_rows, (c + 1) * geometry.num_rows)
                col_patches = [patches[i] for i in col_indices]
                mirrored.extend(col_patches[::-1])  # vertical flip
            patches = mirrored

            for i, p in enumerate(patches, start=1):
                p.index = i

    return patches, white_point_xyz, black_point_xyz



# -------- File Writers --------
class PatchFileWriter:
    """Abstract base class for writers that emit patch measurements to files.

    Subclasses must implement write().

    Attributes provided on init:
      - filename: output path
      - patches: list of PatchInfo
      - column_blocks: list of output blocks (e.g. ['RGB','XYZ','LAB'])
      - white_point_xyz, black_point_xyz: derived XYZ points
      - num_rows, num_cols, row_labels, col_labels: grid metadata
    """

    def __init__(
        self,
        filename: str,
        patches: List[PatchInfo],
        column_blocks: List[str],
        white_point_xyz,
        black_point_xyz,
        num_rows: int,
        num_cols: int,
        row_labels: List[str],
        col_labels: List[str]
    ):
        self.filename = filename
        self.patches = patches
        self.column_blocks = column_blocks
        self.white_point_xyz = white_point_xyz
        self.black_point_xyz = black_point_xyz
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.row_labels = row_labels
        self.col_labels = col_labels

    def write(self):
        raise NotImplementedError("Subclasses must implement write()")


class TI1Writer(PatchFileWriter):
    """Write Argyll .ti1 file.

    TI1 output supports RGB and/or XYZ fields (no Lab). This writer uses
    the provided white_point_xyz/black_point_xyz and the column_blocks to
    format the output. The file follows CTI1 minimal structure.
    """

    def write(self):
        created_ts = make_created_timestamp(time.localtime())
        with open(self.filename, 'w', encoding='utf8') as fh:
            fh.write("CTI1\n\n")
            fh.write('DESCRIPTOR "Argyll Calibration Target chart information 1"\n')
            fh.write('ORIGINATOR "read_image_patch_colors.py (colormath backend)"\n')
            fh.write(f'CREATED "{created_ts}"\n')
            fh.write('ACCURATE_EXPECTED_VALUES "true"\n')

            # Write approximate white point if available
            if self.white_point_xyz is not None:
                Xw, Yw, Zw = [float(v) for v in self.white_point_xyz.tolist()]
                fh.write(f'APPROX_WHITE_POINT "{Xw:.6f} {Yw:.6f} {Zw:.6f}"\n')
            else:
                fh.write('APPROX_WHITE_POINT "0.000000 0.000000 0.000000"\n')

            fh.write('COLOR_REP "iRGB"\n\n')

            # Build headers according to requested column blocks
            headers = ["SAMPLE_ID"]
            for blk in self.column_blocks:
                if blk == 'RGB':
                    headers += ['RGB_R', 'RGB_G', 'RGB_B']
                elif blk == 'XYZ':
                    headers += ['XYZ_X', 'XYZ_Y', 'XYZ_Z']

            fh.write(f'NUMBER_OF_FIELDS {len(headers)}\nBEGIN_DATA_FORMAT\n')
            fh.write(" ".join(headers) + "\nEND_DATA_FORMAT\n\n")

            fh.write(f'NUMBER_OF_SETS {len(self.patches)}\nBEGIN_DATA\n')
            for p in self.patches:
                row_vals = [str(p.index)]
                for blk in self.column_blocks:
                    if blk == 'RGB':
                        row_vals += [f"{v:.6f}" for v in p.rgb_percent]
                    elif blk == 'XYZ':
                        row_vals += [f"{float(v):.6f}" for v in p.xyz100]
                fh.write(" ".join(row_vals) + "\n")
            fh.write('END_DATA\n')

            fh.write("CTI1\n\n")
            fh.write('DESCRIPTOR "Argyll Calibration Target chart information 1"\n')
            fh.write('ORIGINATOR "Argyll targen"\n')
            fh.write('DENSITY_EXTREME_VALUES "8"\n')
            fh.write(f'CREATED "{created_ts}"\n\n')
            fh.write('NUMBER_OF_FIELDS "7"\n')
            fh.write('BEGIN_DATA_FORMAT\n')
            fh.write("INDEX RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z\n")
            fh.write("END_DATA_FORMAT\n\n")
            fh.write('NUMBER_OF_SETS "8"\n')
            fh.write('BEGIN_DATA\n')
            fh.write("0 0.00000 49.99527 43.75000 9.333898 12.12691 15.35996\n")
            fh.write("1 0.00000 50.00123 50.00000 9.282530 11.93912 15.71923\n")
            fh.write("2 0.00000 40.08532 59.90998 9.652263 11.23035 20.85225\n")
            fh.write("3 0.00000 0.00000 78.36374 6.278528 5.962001 12.87821\n")
            fh.write("4 30.55980 50.01120 0.00000 14.44452 18.46006 7.172505\n")
            fh.write("5 0.00000 32.57137 0.00000 5.299583 7.184013 4.303697\n")
            fh.write("6 100.0000 0.00000 0.00000 35.19251 20.43654 4.643321\n")
            fh.write("7 0.00000 0.00000 21.87299 3.801583 3.909116 4.506199\n")
            fh.write('END_DATA\n')

            fh.write("CTI1\n\n")
            fh.write('DESCRIPTOR "Argyll Calibration Target chart information 1"\n')
            fh.write('ORIGINATOR "Argyll targen"\n')
            fh.write('DEVICE_COMBINATION_VALUES "9"\n')
            fh.write(f'CREATED "{created_ts}"\n\n')
            fh.write('NUMBER_OF_FIELDS "7"\n')
            fh.write('BEGIN_DATA_FORMAT\n')
            fh.write("INDEX RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z\n")
            fh.write("END_DATA_FORMAT\n\n")
            fh.write('NUMBER_OF_SETS "9"\n')
            fh.write('BEGIN_DATA\n')
            fh.write("0 100.0000 100.0000 100.0000 92.38129 96.32111 81.69708\n")
            fh.write("1 0.00000 100.0000 100.0000 29.15368 33.35974 49.04946\n")
            fh.write("2 100.0000 0.00000 100.0000 36.32309 22.82926 20.84792\n")
            fh.write("3 0.00000 0.00000 100.0000 6.999917 6.572461 14.97645\n")
            fh.write("4 100.0000 100.0000 0.00000 75.23812 85.69790 15.35314\n")
            fh.write("5 0.00000 100.0000 0.00000 16.87349 25.12038 10.16533\n")
            fh.write("6 100.0000 0.00000 0.00000 35.19251 20.43654 4.643321\n")
            fh.write("7 0.00000 0.00000 0.00000 4.312931 5.032890 3.693281\n")
            fh.write("8 50.00000 50.00000 50.00000 24.83442 26.40655 22.95545\n")
            fh.write('END_DATA\n')

class TI2Writer(PatchFileWriter):
    """Write Argyll .ti2 file.

    TI2 includes sample ID, sample location, and the chosen color columns
    (RGB, XYZ, LAB). It also emits metadata such as steps/passes and index patterns.
    """

    def write(self):
        created_ts = make_created_timestamp(time.localtime())

        # Helper to detect label type (alphabetic vs numeric)
        def is_alpha_labels(lbls: List[str]) -> bool:
            return bool(re.fullmatch(r"[A-Z]+", lbls[0]))

        if is_alpha_labels(self.row_labels):
            strip_index_pattern = 'A-Z, A-Z'
        else:
            strip_index_pattern = '0-9,@-9;1-99'

        if is_alpha_labels(self.col_labels):
            patch_index_pattern = 'A-Z, A-Z'
        else:
            patch_index_pattern = '0-9,@-9;1-99'

        with open(self.filename, 'w', encoding='utf8') as fh:
            fh.write("CTI2\n\n")
            fh.write('DESCRIPTOR "Argyll Calibration Target chart information 2"\n')
            fh.write('ORIGINATOR "read_image_patch_colors.py (colormath backend)"\n')
            fh.write(f'CREATED "{created_ts}"\n')
            fh.write('ACCURATE_EXPECTED_VALUES "true"\n')

            # Write approximate white point if available
            if self.white_point_xyz is not None:
                Xw, Yw, Zw = [float(v) for v in self.white_point_xyz.tolist()]
                fh.write(f'APPROX_WHITE_POINT "{Xw:.6f} {Yw:.6f} {Zw:.6f}"\n')
            else:
                fh.write('APPROX_WHITE_POINT "0.000000 0.000000 0.000000"\n')

            fh.write('COLOR_REP "iRGB"\n')

            # Number of rows/cols and index patterns
            fh.write(f'STEPS_IN_PASS "{self.num_rows}"\n')
            fh.write(f'PASSES_IN_STRIPS2 "{self.num_cols}"\n')
            fh.write(f'STRIP_INDEX_PATTERN "{strip_index_pattern}"\n')
            fh.write(f'PATCH_INDEX_PATTERN "{patch_index_pattern}"\n')
            fh.write('INDEX_ORDER "STRIP_THEN_PATCH"\n\n')

            # Headers include SAMPLE_ID and SAMPLE_LOC
            headers = ["SAMPLE_ID", "SAMPLE_LOC"]
            for blk in self.column_blocks:
                if blk == 'RGB':
                    headers += ['RGB_R', 'RGB_G', 'RGB_B']
                elif blk == 'XYZ':
                    headers += ['XYZ_X', 'XYZ_Y', 'XYZ_Z']
                elif blk == 'LAB':
                    headers += ['LAB_L', 'LAB_A', 'LAB_B']

            fh.write(f'NUMBER_OF_FIELDS {len(headers)}\nBEGIN_DATA_FORMAT\n')
            fh.write(" ".join(headers) + "\nEND_DATA_FORMAT\n\n")

            fh.write(f'NUMBER_OF_SETS {len(self.patches)}\nBEGIN_DATA\n')
            for p in self.patches:
                row_vals = [str(p.index), f'"{p.label}"']
                for blk in self.column_blocks:
                    if blk == 'RGB':
                        row_vals += [f"{v:.6f}" for v in p.rgb_percent]
                    elif blk == 'XYZ':
                        row_vals += [f"{float(v):.6f}" for v in p.xyz100]
                    elif blk == 'LAB':
                        row_vals += [f"{float(v):.6f}" for v in p.lab]
                fh.write(" ".join(row_vals) + "\n")
            fh.write('END_DATA\n')


class CSVWriter(PatchFileWriter):
    """Write a space-separated CSV file containing the same columns as TI2.

    This writer produces a plain text file with the header row followed by data rows.
    """

    def write(self):
        headers = ["SAMPLE_ID", "SAMPLE_LOC"]
        for blk in self.column_blocks:
            if blk == 'RGB':
                headers += ['RGB_R', 'RGB_G', 'RGB_B']
            elif blk == 'XYZ':
                headers += ['XYZ_X', 'XYZ_Y', 'XYZ_Z']
            elif blk == 'LAB':
                headers += ['LAB_L', 'LAB_A', 'LAB_B']

        with open(self.filename, 'w', encoding='utf8') as fh:
            fh.write(" ".join(headers) + "\n")
            for p in self.patches:
                row_vals = [str(p.index), p.label]
                for blk in self.column_blocks:
                    if blk == 'RGB':
                        row_vals += [f"{v:.6f}" for v in p.rgb_percent]
                    elif blk == 'XYZ':
                        row_vals += [f"{float(v):.6f}" for v in p.xyz100]
                    elif blk == 'LAB':
                        row_vals += [f"{float(v):.6f}" for v in p.lab]
                fh.write(" ".join(row_vals) + "\n")


# -------- Main Entry Point --------
# Process the image and write .ti1, .ti2 and .csv output files.
# Input: args (argparse.Namespace) with all CLI flags.
# Behavior: validates labels, computes geometry, samples patches, converts colours, writes outputs.
def process_image_to_files(args):
    # Parse and validate label ranges
    row_labels, row_pad = parse_label_range(args.row_labels)
    col_labels, col_pad = parse_label_range(args.col_labels)
    if len(row_labels) != args.num_rows or len(col_labels) != args.num_cols:
        raise ValueError("Label counts do not match num_rows/num_cols")

    print("Row/column label counts match grid dimensions. Continuing...\n")

    # Parse output color space tokens and validate
    included = [tok.strip().upper() for tok in args.output_color_space.split(',') if tok.strip()]
    valid = {"RGB", "XYZ", "LAB"}
    for tok in included:
        if tok not in valid:
            raise ValueError(
                f"Invalid output_color_space token '{tok}'. At least one of RGB, XYZ, LAB is required"
            )

    # Load image and compute geometry
    img16, _ = load_image_as_16bit_rgb(args.image)
    image_color_space = args.image_color_space.lower() if args.image_color_space else 'srgb'
    geometry = compute_grid_geometry(args)

    # Extract measurements
    patches, white_point_xyz, black_point_xyz = extract_patch_data(
        img16,
        geometry,
        row_labels,
        col_labels,
        row_pad,
        col_pad,
        args.patch_label_order,
        image_color_space,
        args.sample_mode,
        args.output_order,
        args.mirror_output,
        args.rotate_grid
    )

    base_out = os.path.splitext(os.path.basename(args.image))[0]

    ti1_file = args.output if args.output else f"{base_out}.ti1"
    ti2_file = args.output if args.output else f"{base_out}.ti2"
    csv_file = f"{base_out}.csv"

    # Instantiate and write each requested file type using common writer interface
    for writer_cls, fname in [
        (TI1Writer, ti1_file),
        (TI2Writer, ti2_file),
        (CSVWriter, csv_file),
    ]:
        writer = writer_cls(
            fname,
            patches,
            included,
            white_point_xyz,
            black_point_xyz,
            args.num_rows,
            args.num_cols,
            row_labels,
            col_labels,
        )
        writer.write()
        print("Wrote:", fname)

    print("Done.")


# ---------- CLI ----------
# Build argparse parser with expected command line options.
def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Read patch colours from grid image with label metadata.")
    p.add_argument('--image', '-i', required=True, help='Image file path')
    p.add_argument('--image_color_space', default='srgb', choices=['srgb', 'adobergb'],
                   help='Input image device colour space (default: srgb)')
    p.add_argument('--patch_first_xy', required=True, type=parse_xy,
                   help='X,Y of first patch (top-left)')
    p.add_argument('--patch_last_xy', required=True, type=parse_xy,
                   help='X,Y of last patch (bottom-right)')
    p.add_argument('--patch_width_height_ratio', required=True, type=float,
                   help='Width / Height ratio for a single patch')
    p.add_argument('--num_cols', required=True, type=int, help='Number of columns in grid')
    p.add_argument('--num_rows', required=True, type=int, help='Number of rows in grid')
    p.add_argument('--sample_fraction', type=float, default=DEFAULT_SAMPLE_FRACTION,
                   help='Fraction (0<f≤0.6) of smaller patch side to sample')
    p.add_argument('--row_labels', required=True, help='Range for rows, e.g. "1-20" or "A-D"')
    p.add_argument('--col_labels', required=True, help='Range for columns, e.g. "1-10" or "A-C"')
    p.add_argument('--patch_label_order', required=True, choices=["col_then_row", "row_then_col"])
    p.add_argument('--output_color_space', required=True,
                   help='Comma-separated sequence of color spaces to include: RGB, XYZ, LAB')
    p.add_argument(
        '--sample_mode',
        default='mad',
        choices=['mean', 'median', 'mad'],
        help=(
            "Sampling aggregation method:\n"
            "  mean   – plain arithmetic mean (sensitive to outliers)\n"
            "  median – pure median (very robust, but may bias bright/dark values)\n"
            "  mad    – robust mean using MAD-based sigma clipping (default)"
        )
    )
    p.add_argument('--output_order', default='row_major', choices=['row_major', 'column_major'], help='Output patch data one row at a time (row_major), or one column at a time (column_major). Applied before rotate_grid and mirror_output')
    p.add_argument('--mirror_output', action='store_true', help='Output patch data with reversed traversal. If row_major then reversed column traversal. If column_major then reversed row traversal. Applied after output_order and rotate_grid.')
    p.add_argument(
        '--rotate_grid',
        default=0,
        type=int,
        choices=[0, 90, 180, 270],
        help='Rotate the patch grid by 0 (default), 90, 180, or 270 degrees clockwise before output. Applied after output_order and before mirror_output.'
    )
    p.add_argument('--output', help='Optional base name for output TI1/TI2/CSV files')
    p.add_argument('--debug', action='store_true',
                   help='Enable diagnostic debug printing for the first patch')
    return p


def main():
    global DEBUG, DEBUG_AVG_COUNTER, DEBUG_SAMPLE_COUNTER

    parser = build_arg_parser()
    args = parser.parse_args()

    DEBUG = args.debug
    DEBUG_AVG_COUNTER = 0
    DEBUG_SAMPLE_COUNTER = 0
    try:
        process_image_to_files(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write("Fatal error: " + str(e) + "\n")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
