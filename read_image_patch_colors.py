#!/usr/bin/env python3
"""
read_image_patch_colors.py
Version: 1.5.0

================================================================================
READ_IMAGE_PATCH_COLORS.PY — DOCUMENTATION
================================================================================

OVERVIEW
--------
This program extracts color values from a rectangular grid of color patches in
an image, computes colorimetric values (RGB percentages, XYZ, Lab and approximate
white point), applies row/column labeling rules, and writes three output files:

    • ArgyllCMS .ti1 file
    • ArgyllCMS .ti2 file
    • CSV file (space-separated)

It supports sampling in mean, median or MAD-based sigma clipping mode, configurable
patch geometry, numeric or alphabetic labeling patterns, and RGB/XYZ/Lab output
combinations.

The script is designed for use in color-management workflows where printed color
targets are scanned or photographed and converted into Argyll measurement files.
Alternatively, convert original calibration target images to attain `.ti1`, `.ti2`
and `.csv` files with their patch color data.

Example use cases, using the image of a reference target:

1. Create a `.ti1` file so that one may use ArgyllCMS `printtarg` command and
   generate a target using the colors of the image, which then can be used with
   ArgyllCMS `chartread` and `colprof` to create a printer profile.
    - Workflow:
      Read image with script → input `.ti1` to `printtarg` → Print generated Target
      image without color management → input `.ti2` to `chartread` then read target
      paches with colorimeter → input `.ti3` to `colprof` → output profile.

2. Create a `.ti2` file so that one may use ArgyllCMS `chartread`, and use a
   colorimeter to read the color values using the original SpyderPrint Targets,
   which then can be used with ArgyllCMS `colprof` to create a printer profile.
    - Workflow:
    Read image with script → Print original SpyderPrint Target image without color
    management → input `.ti2` to `chartread` then read target paches with colorimeter
    → input `.ti3` to `colprof` → output profile.

3. Create a `.ti1` file so that one may use ArgyllCMS `fakeread` command using
   colors of the image, which then can be used with ArgyllCMS `colprof` to create
   a simulated profile.
    - Workflow:
    Read image with script → input `.ti1` to `fakeread` → input `.ti3` to `colprof`
    → output simulated profile.

--------------------------------------------------------------------------------
IMAGE INPUT
-----------
Accepted image types: any PIL-compatible raster file
Color depth: handled as 16-bit RGB internally


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
    • "mad" (default)
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

Valid range:  0.0 - 100.0

Color conversions (XYZ, Lab) are computed through ArgyllCMS xicclu command with
selected color icc profile (default `sRGB.icm`).
Absolute colorimetric intent with D50 illuminant is used, same as how ArgyllCMS
targen creates patch colors, and is expected by printtarg and chartread.
Targen uses a sRGB-like model for generating XYZ and LAB D50, unless another
profile i provided. This is why `sRGB.icm` color space profile is used by default,
but may be overridden by choosing a different profile with `--pre_cond_profile`.

XYZ output:

- CIEXYZ (CIE 1931 2° standard observer)
- D50 white reference (Profile Connection Space in ICC)
- Absolute colorimetric intent
- White scaled to Y=100, unbounded positive, values 0-100
- Numeric precision: six decimal places.

Lab output:

- CIE L*a*b* (1976), computed using **D50** reference white (Xn, Yn, Zn)
- Absolute colorimetric intent
- L in range ~0..100, a and b range typical -128..128.
- Numeric precision: six decimal places.

- All colorimetric output is **D50** and **linear** (no gamma applied in XYZ or Lab).

The patch with the highest Y (XYZ_Y) is stored as APPROX_WHITE_POINT. Note that XYZ
Absolute colorimetric D50 is almost identical to XYZ Relative colorimetric D65,
which can confuse some users to think that the APPROX_WHITE_POINT is computed
using D65. This is not the case.

--------------------------------------------------------------------------------
OUTPUT FILES
------------
Three files are produced:

1. TI1 FILE (.ti1)
------------------
Contains:
    • Header (descriptor, originator, timestamp)
    • COLOR_REP (always "iRGB")
    • ACCURATE_EXPECTED_VALUES flag
    • Data fields (SAMPLE_ID, RGB, optionally XYZ)
    • Patch data in Argyll .ti1 format (SAMPLE_ID, RGB and optionally XYZ)
    • Dynamically generated:
        • APPROX_WHITE_POINT
        • WHITE_COLOR_PATCHES
        • BLACK_COLOR_PATCHES
        • COMP_GREY_STEPS
        • SINGLE_DIM_STEPS
        • NUMBER_OF_FIELDS
        • NUMBER_OF_SETS
        • DENSITY_EXTREME_VALUES table
        • DEVICE_COMBINATION_VALUES table

Tags containing a dynamically calculated number are omitted if 0 (zero).

COMP_GREY_STEPS
Number represents count of patches with equal RGB channel values.
RGB Balance Check:
    - RGB Tolerance: 0.1%. Very tight tolerance requiring near-perfect R=G=B balance
    - Lightness Range: 1.0% to 99.0%. Excludes pure black and white patches.
    - Logic: A patch is neutral if all RGB channels are balanced within 0.1% and
      fall within the specified lightness range.

SINGLE_DIM_STEPS
This number is estimated by reading single-dimension color ramps where two RGB
channels remain static while the third creates a ramp pattern. The number
represents the most commonly detected number of steps for single-channel ramps.
A single-channel ramp is defined as a sequence where:
    - Two RGB channels remain constant (within tolerance)
    - One RGB channel changes progressively (positive or negative direction)
    - At least 3 patches are required to form a valid ramp
    - Values are approximately evenly spaced (uniform spacing between steps)
    - Patches can be in any order (randomized)

2. TI2 FILE (.ti2)
------------------
Contains:
    • Header (descriptor, originator, timestamp)
    • STRIP and PATCH index patterns inferred from row/col labels
    • STEPS_IN_PASS     = num_rows
    • PASSES_IN_STRIPS2 = num_cols
    • Index order       = STRIP_THEN_PATCH
    • Dynamically generated:
        • APPROX_WHITE_POINT
        • NUMBER_OF_FIELDS
        • NUMBER_OF_SETS
    • Complete patch data (SAMPLE_ID, SAMPLE_LOC, RGB, XYZ, Lab)

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
    • ArgyllCMS

If missing, the script prints installation instructions and exits.

--------------------------------------------------------------------------------
ERROR CONDITIONS
----------------
The script stops with a readable message when:
    • Coordinates and Row/Column numbers cause measurement areas to be missaligned
      with patch centres of the target image.
    • Labels do not match row/column counts
    • Output color space tokens are invalid
    • sample_fraction exceeds allowed range
    • Grid geometry cannot be computed
    • Any I/O error occurs during file writing
    • File path to image, pre_cond_profile, or output directory is invalid

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
                                Input image must be an image target made for printing.
--pre_cond_profile              ICC/ICM profile file path, used as device
                                pre-conditioning profile. Default: 'sRGB.icm'.
                                This adapts/shapes the XYZ and LAB values to a known
                                device/printer profile or color space profile.
                                This also affects the approx. white point value.
--patch_first_xy                "X,Y" coordinates (floats or ints) of the ***centre***
                                of the first patch (top-left patch). Origin is top-left
                                of the image.
--patch_last_xy                 "X,Y" coordinates (floats or ints) of the ***centre***
                                of the last patch (bottom-right patch). Accepts integers
                                or floats.
                                Note: In some cases some pathces on last row of image
                                may not exist. In these cases to ccordinates to the
                                last row and column, as if it would exist, must be
                                provided. After generating of .ti1, .ti2 and .csv
                                files the last patches must be manually removed, if
                                output must be exact like original image.
                                If patch_label_order is col_then_row then the patces
                                to remove are in sequence at the end of the list, but
                                if row_then_col is used, then you would need to search
                                for the patch lables (SAMPLE_LOC) to locate them and
                                remove them.
--patch_width_height_ratio      Width ÷ height (W/H) ratio of a single patch.
                                    Examples:
                                      1.0   → square patches
                                      1.4   → width is 1.4× height
--num_cols                      Number of columns in the grid.
--num_rows                      Number of rows in the grid.
--row_labels
--col_labels
                                Define the label sequences. Must match num_rows / num_cols
                                exactly.
                                Supported formats:
                                • Numeric: 1-27, 03-15, 0001-0120 (zero-padding preserved)
                                • Alphabetic: A-Z, A-AA, BQ-CF, up to ZZ
                                Notes:
                                • Alphabetic and numeric are mutually exclusive between
                                row and column: If rows use alphabetic labels, columns
                                must be numeric, and vice-versa.
--patch_label_order             Determines how the patch label is formed:
                                  col_then_row → "<col><row>"
                                  row_then_col → "<row><col>"
--output_color_space           Comma-separated list. Allowed tokens: RGB, XYZ, LAB
                                Specifies **which colour spaces** to output in generated
                                file, and in what order.
                                TI1 only includes RGB and/or XYZ. TI2/CSV includes any
                                defined token.
                                Example tokens:
                                  RGB,XYZ       →   TI1 includes RGB,XYZ; TI2/CSV includes RGB,XYZ
                                  XYZ,LAB,RGB   →   TI1 includes XYZ,RGB; TI2/CSV includes XYZ,LAB,RGB
                                Output columns follow the specified sequence.
--sample_fraction               Fraction of patch to sample; default 0.50 (50%).
                                The sampling square is centred on the patch centre and
                                is clamped to at least 3×3 pixels and at most 60% of patch size.
--sample_mode                    "mean" (default) or "median" for robust sampling
--output                         [Optional] output filename (defaults to <imagebasename>.ti1/ti2/csv)

================================================================================
COLOR / SCALING CONVENTIONS IN OUTPUT FILES (Argyll-compatible)
================================================================================
- RGB are written as **iRGB** percentage values: 0.000000 .. 100.000000
  (Argyll convention; device nominal percentages).

- CIEXYZ in output:
  * CIEXYZ (CIE 1931 2° standard observer)
  * **D50** white reference (Profile Connection Space in ICC)
  * **Scaled so that white has Y = 100.000000**
  * Numeric precision: six decimal places.

- Lab in output:
  * CIE L*a*b* (1976), computed using **D50** reference white (Xn, Yn, Zn)
  * L in range ~0..100, a and b around typical -128..128 ranges.
  * Numeric precision: six decimal places.

- All colorimetric output is **D50** and **linear** (no gamma applied in XYZ or Lab).
  The conversion pipeline in this script applies conversions according to
  the selected `--pre_cond_profile` using ArgyllCMS xicclu.

================================================================================
NOTES ON ACCURACY
================================================================================
- This script seeks to minimize rounding/quantization error:
  * If the image has 16-bit/channel data it will be used directly.
  * If the image is 8-bit/channel (most typical images), it will be scaled up
    exactly to 16-bit (value * 257) before colour conversion to reduce rounding
    / mapping error in conversions.
  * Conversion from RGB -> XYZ -> Lab is done using ArgyllCMS xicclu with
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
4. The sampled RGB values are averaged (MAD-based sigma) in 16-bit integer space
   (8-bit inputs are scaled to 16-bit via *257).
5. Colour conversions are applied using xicclu:
   - RGB -> CIEXYZ XYZ (D50)
   - RGB -> Lab (D50)

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
      --patch_first_xy 72,51 \
      --patch_last_xy 1250,934 \
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
import subprocess
import shutil
import importlib.util
from dataclasses import dataclass
from typing import List, Tuple, Sequence, Optional

# ---------------- Preflight ----------------
# Check that required Python packages are installed.
# Input: None
# Output: list of missing package names (list[str]).
# Packages checked: numpy, Pillow
def check_full_packages() -> list[str]:
    missing = []

    # 1. Check ArgyllCMS by checking whether xicclu exists in PATH
    if shutil.which("xicclu") is None:
        missing.append("ArgyllCMS (xicclu not found in PATH)")

    # 2. Check Python packages
    if importlib.util.find_spec("numpy") is None:
        missing.append("numpy")

    # Pillow must be checked via PIL.Image, not PIL
    if importlib.util.find_spec("PIL.Image") is None:
        missing.append("Pillow")

    return missing


# Print platform-specific pip install instructions for missing packages.
def print_install_instructions(packages: list[str]) -> None:
    system = platform.system()
    print("\nMissing components detected:\n")
    for pkg in packages:
        print(" -", pkg)

    print("\nTo install:")

    if any("ArgyllCMS" in p for p in packages):
        if system == "Darwin":
            print("\nInstall ArgyllCMS on macOS:")
            print("  brew install argyll-cms")
        elif system == "Windows":
            print("\nInstall ArgyllCMS on Windows:")
            print("  Download installer:")
            print("  https://www.argyllcms.com/")
        else:
            print("\nInstall ArgyllCMS for Linux:")
            print("  sudo apt install argyll")

    # Python packages
    py_pkgs = [p for p in packages if "ArgyllCMS" not in p]
    if py_pkgs:
        if system == "Windows":
            print(f"\nPython packages:")
            print(f"  python -m pip install {' '.join(py_pkgs)}")
        else:
            print(f"\nPython packages:")
            print(f"  python3 -m pip install {' '.join(py_pkgs)}")


missing = check_full_packages()
if missing:
    print_install_instructions(missing)
    sys.exit(1)


# ---------------- Imports AFTER preflight ----------------
from PIL import Image
import numpy as np

# ---------- Constants ----------
DEBUG = False
DEBUG_PRINT_LIMIT = 10
DEBUG_AVG_COUNTER = 0
EXTREME_OUTLIER_PATCHES = 0
EXTREME_OUTLIER_LIMIT = 5
DEBUG_SAMPLE_COUNTER = 0
DEFAULT_SAMPLE_FRACTION = 0.50
EPS = (6.0 / 29.0) ** 3   # Threshold for linearization in Lab
K = 24389.0 / 27.0        # Linear coefficient in Lab conversion

# Match ANY float (scientific, signed, etc.)
FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

# Match 3 consecutive floats
TRIPLE_RE = re.compile(rf"({FLOAT_RE}\s+{FLOAT_RE}\s+{FLOAT_RE})")

# ------------------------------------------------
# RGB lists (0..100, as floats)
# ------------------------------------------------
density_extr_val_list = [
    (0.00000, 49.99527, 43.75000),
    (0.00000, 50.00123, 50.00000),
    (0.00000, 40.08532, 59.90998),
    (0.00000, 0.00000, 78.36374),
    (30.55980, 50.01120, 0.00000),
    (0.00000, 32.57137, 0.00000),
    (100.0000, 0.00000, 0.00000),
    (0.00000, 0.00000, 21.87299),
]

device_combi_val_list = [
    (100.0000, 100.0000, 100.0000),
    (0.00000, 100.0000, 100.0000),
    (100.0000, 0.00000, 100.0000),
    (0.00000, 0.00000, 100.0000),
    (100.0000, 100.0000, 0.00000),
    (0.00000, 100.0000, 0.00000),
    (100.0000, 0.00000, 0.00000),
    (0.00000, 0.00000, 0.00000),
    (50.00000, 50.00000, 50.00000),
]

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


# Load an image and return a 16-bit RGB numpy array (HxWx3, uint16).
# Input: img_path (str), assumes images is RGB and allows no other mode.
#        Reads RGB 3x8-bit pixels, or 3x16 pixels, true color.
# Output: tuple (arr16: np.ndarray HxWx3 uint16, bit_depth: int 8|16)
def load_image_as_16bit_rgb(img_path: str) -> Tuple[np.ndarray, int]:
    try:
        im = Image.open(img_path)
    except Exception as e:
        print("Error: Failed to open image:", e)
        sys.exit(1)

    # Only allow RGB mode
    if im.mode != "RGB":
        print(f"Error: Unsupported image mode '{im.mode}'. Only RGB images are allowed.")
        sys.exit(1)
    # Image mode
    mode = im.mode

    #debug print
    if DEBUG:
        print(f"Image has mode: '{im.mode}'. Only RGB images are suppored.")

    # Convert to RGB (ensures no alpha, no grayscale)
    im = im.convert("RGB")
    arr = np.array(im)

    if arr.dtype == np.uint8:
        # Expand 8-bit to full 16-bit range (value * 257)
        arr16 = (arr.astype(np.uint32) * 257).astype(np.uint16)
        return arr16, 8

    elif arr.dtype == np.uint16:
        # Already 16-bit per channel
        return arr.astype(np.uint16), 16

    else:
        # Floating point or other → normalize 0..1 then scale
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

    # ============================================================
    # 1. Compute output value depending on mode (but DO NOT RETURN)
    # ============================================================
    # ----- simple mean -----
    if mode == 'mean':
        result = arr.mean(axis=0)

    # ----- simple median -----
    elif mode == 'median':
        result = np.median(arr, axis=0)

    else:  # mode == 'mad'
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
        result = np.array([r, g, b], dtype=np.float64)

        # Debug native and robust mean for first 10 MAD patches
        if DEBUG and mode == 'mad' and DEBUG_AVG_COUNTER < DEBUG_PRINT_LIMIT:
            patch_num = DEBUG_AVG_COUNTER + 1
            print(f"--- PATCH {patch_num} ---")
            print(f"Naive mean (16-bit): {arr.mean(axis=0)}")
            print(f"Robust mean (16-bit): {result}")
            # Show MAD / sigma info per channel
            for i, ch in enumerate(["R", "G", "B"]):
                med = np.median(arr[:, i])
                mad = np.median(np.abs(arr[:, i] - med))
                sigma = max(1.4826 * mad, 1.0)
                mask = np.abs(arr[:, i] - med) <= 3 * sigma
                rejected = (~mask).sum()
                print(f"{ch}-channel: median={med:.1f}, MAD={mad:.1f}, sigma={sigma:.1f}, outliers removed={rejected}")
            DEBUG_AVG_COUNTER += 1

    # ------------------------------------------------------------
    #  EXTREME OUTLIER DETECTION (10% deviation rule)
    #  We compare every pixel value with the chosen aggregate:
    #     mad → result
    # ------------------------------------------------------------
    global EXTREME_OUTLIER_PATCHES

    # Compute MAD-based reference (per channel) for *all* modes
    def mad_ref(values: np.ndarray, min_sigma: float = 1.0):
        med = np.median(values)
        mad = np.median(np.abs(values - med))
        sigma = max(1.4826 * mad, min_sigma)
        return med, sigma

    # Per-channel median + sigma
    med_r, sig_r = mad_ref(arr[:, 0])
    med_g, sig_g = mad_ref(arr[:, 1])
    med_b, sig_b = mad_ref(arr[:, 2])

    med_vec = np.array([med_r, med_g, med_b], dtype=np.float64)
    sigma_vec = np.array([sig_r, sig_g, sig_b], dtype=np.float64)

    # Define "extreme" as 10% of median or 5 units minimum
    abs_thresh = np.maximum(np.abs(med_vec) * 0.10, 5.0)

    # Pixels where ANY channel deviates too far from the MAD-median
    diffs = np.abs(arr - med_vec) > abs_thresh
    has_extreme = np.any(np.any(diffs, axis=1))

    if has_extreme:
        EXTREME_OUTLIER_PATCHES += 1

    return result


# Convert 16-bit RGB to Argyll iRGB percentage (0.0 .. 100.0)
# Input: rgb16 (3 values 0..65535)
# Output: tuple of floats in 0..100
def rgb16_to_argyll_percent(rgb16: Sequence[float]) -> Tuple[float, float, float]:
    factor = 100.0 / 65535.0
    return float(rgb16[0] * factor), float(rgb16[1] * factor), float(rgb16[2] * factor)


def parse_xicclu_output(text: str):
    """
    Bulletproof parser for xicclu output.

    Strategy:
      • For each non-empty line:
          - Find *all* float-triplets in the line
          - Use the LAST triplet (this is always XYZ or LAB)
      • Return list of 3-tuples of floats

    If debug=True, prints how each line is parsed.
    """

    results = []

    if DEBUG:
        print(f"\n\nConfirming filtering of output from xicclu in parse_xicclu_output:")

    for line in text.splitlines():
        raw_line = line.rstrip("\n")
        line = raw_line.strip()

        if not line:
            continue

        # Find ALL 3-float sequences
        triples = TRIPLE_RE.findall(line)

        if not triples:
            if DEBUG:
                print(f"[SKIP] No triples found in line: '{raw_line}'")
            continue

        # Use LAST triple
        trip = triples[-1]

        parts = trip.split()
        if len(parts) != 3:
            if DEBUG:
                print(f"[SKIP] Invalid triple: '{trip}' in line: '{raw_line}'")
            continue

        tup = tuple(map(float, parts))
        results.append(tup)

        if DEBUG:
            print(f"[OK] {tup}  ← from: '{raw_line}'")

    return results


def convert_color_batch(rgb_percent_list, icc_profile_path):
    """
    Batch-convert a list of RGB percentages into XYZ and Lab using xicclu,
    with only one subprocess call for XYZ and one for Lab.

    rgb_percent_list: list of (R,G,B) floats 0..100
    Returns:
        xyz_list: list of (X,Y,Z)
        lab_list: list of (L,a,b)
    """

    # ------------------------------------------------
    # Build xicclu input block
    # ------------------------------------------------
    xicclu_input = "\n".join(
        f"{rgb[0]} {rgb[1]} {rgb[2]}" for rgb in rgb_percent_list
    ) + "\n\n"

    # ------------------------------------------------
    # Run xicclu for XYZ (use -i a = absolute colorimetric intent with D50, same as targen does)
    # ------------------------------------------------
    cmd_xyz = ["xicclu", "-ff", "-ia", "-s100", "-pX", icc_profile_path]
    out_xyz = subprocess.run(
        cmd_xyz, input=xicclu_input, capture_output=True, text=True
    ).stdout

    # ------------------------------------------------
    # Run xicclu for Lab (use -i a = absolute colorimetric intent with D50, same as targen does)
    # ------------------------------------------------
    cmd_lab = ["xicclu", "-ff", "-ia", "-s100", "-pl", icc_profile_path]
    out_lab = subprocess.run(
        cmd_lab, input=xicclu_input, capture_output=True, text=True
    ).stdout

    # Parse XYZ
    xyz_list = parse_xicclu_output(out_xyz)
    # Parse Lab
    lab_list = parse_xicclu_output(out_lab)

    if len(xyz_list) != len(rgb_percent_list) or len(lab_list) != len(rgb_percent_list):
        raise ValueError("Parsed xicclu batch results do not match input count.")

    return xyz_list, lab_list


def build_patch_records(rgb_percent_list, icc_profile_path):
    """
    Convenience: return a list of per-patch records:
        (index, (R,G,B), (X,Y,Z), (L,a,b))
    """
    xyz_list, lab_list = convert_color_batch(rgb_percent_list, icc_profile_path)

    records = []
    for idx, (rgb, xyz, lab) in enumerate(zip(rgb_percent_list, xyz_list, lab_list)):
        records.append((idx, rgb, xyz, lab))
    return records


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
#   pre_cond_profile : string, color space icc profile file patch
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
    pre_cond_profile: str,
    sample_mode: str,
    output_order: str,
    mirror_output: bool,
    rotate_grid: int
) -> Tuple[List[PatchInfo], np.ndarray, np.ndarray]:

    patches: List[PatchInfo] = []
    patch_idx = 0
    white_point_xyz = None
    black_point_xyz = None
    total_patches = geometry.num_rows * geometry.num_cols
    rgb_percent_list = []
    # mapping from patch label to its RGB index in rgb_percent_list
    rgb_index_by_label = {}

    # ------------------------------------------------------
    # ROW-MAJOR ORDER
    # ------------------------------------------------------
    if output_order == "row_major":

        for row in range(geometry.num_rows):
            for col in range(geometry.num_cols):

                patch_idx += 1
                # -------------------------------
                # Show progress
                # -------------------------------
                percent = (patch_idx / total_patches) * 100
                sys.stdout.write(f"\rProcessing patch {patch_idx}/{total_patches} ({percent:.1f}%)")
                sys.stdout.flush()

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
                rgb_percent_list.append(rgb_percent)
                # Track RGB index by label for later
                rgb_index_by_label[patch_label] = len(rgb_percent_list) - 1

                # ------------------------------------------------------
                # EXCESSIVE OUTLIER SAFETY CHECK
                # ------------------------------------------------------
                if EXTREME_OUTLIER_PATCHES > EXTREME_OUTLIER_LIMIT:
                    print("\n\nError: multiple patch measurements has detected extreme value outliers. "
                          "This may be due to the measurement area not being placed properly within the "
                          "patch areas of the target image. Typically, wrong input in any of the following "
                          "arguments may cause this: '--patch_first_xy', '--patch_last_xy', "
                          "'--patch_width_height_ratio', '--num_cols' or '--num_rows'.\n\n"
                          "Please verify command inputs and try again.\n")
                    sys.exit(1)

                patches.append(
                    PatchInfo(
                        index=patch_idx,
                        label=patch_label,
                        rgb16=avg_rgb16,
                        rgb_percent=rgb_percent,
                        xyz100=(None, None, None),
                        lab=(None, None, None)
                    )
                )

        # ----- Rotate grid before mirroring -----
        if rotate_grid != 0:
            # After building patches, reshape into 2D grid
            patch_grid = np.array(patches, dtype=object).reshape((geometry.num_rows, geometry.num_cols))
            # Apply rotation
            k = rotate_grid // 90  # np.rot90 rotates clockwise
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

                # -------------------------------
                # Show progress
                # -------------------------------
                percent = (patch_idx / total_patches) * 100
                sys.stdout.write(f"\rProcessing patch {patch_idx}/{total_patches} ({percent:.1f}%)")
                sys.stdout.flush()

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
                rgb_percent_list.append(rgb_percent)
                # Track RGB index by label for later
                rgb_index_by_label[patch_label] = len(rgb_percent_list) - 1

                patches.append(
                    PatchInfo(
                        index=patch_idx,
                        label=patch_label,
                        rgb16=avg_rgb16,
                        rgb_percent=rgb_percent,
                        xyz100=(None, None, None),
                        lab=(None, None, None)
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

    # ------------------------------------------------------
    # FINAL STEP: batch-convert all RGB% values using xicclu
    # ------------------------------------------------------
    xyz_list, lab_list = convert_color_batch(rgb_percent_list, pre_cond_profile)

    if len(xyz_list) != len(patches):
        raise RuntimeError("XYZ/Lab list size mismatch with patches.")

    # Assign results to each patch (preserves order, label-based assignment)
    for p in patches:
        idx = rgb_index_by_label[p.label]   # get the original RGB index for this label
        p.xyz100 = xyz_list[idx]
        p.lab    = lab_list[idx]

    white_point_xyz = None
    black_point_xyz = None

    for p in patches:
        Y = p.xyz100[1]
        if white_point_xyz is None or Y > white_point_xyz[1]:
            white_point_xyz = p.xyz100
        if black_point_xyz is None or Y < black_point_xyz[1]:
            black_point_xyz = p.xyz100

    # ------------------------------------------------------
    # DEBUG: Print XYZ/LAB for first 10 patches
    # ------------------------------------------------------
    if DEBUG:
        print("\n\n--- DEBUG: First 10 patches with XYZ/LAB ---")
        for i, p in enumerate(patches[:10], 1):
            print(f"\nPatch {i}: {p.label}")
            print(f"RGB16    : {p.rgb16}")
            print(f"RGB %    : {p.rgb_percent}")
            print(f"XYZ100   : {p.xyz100}")
            print(f"LAB      : {p.lab}")

    print("\nAll patches processed.")

    return patches, white_point_xyz, black_point_xyz


def resolve_icc_profile_path(profile_arg: str, arg_name: str) -> str:
    """
    Resolve and validate ICC/ICM profile path.

    Behavior:
    - If the user did NOT specify the argument at all → return sRGB.icm (use fallback).
    - If the user specified the flag but DID NOT provide a value → python argument error.
    - If they provided a value → validate it.
    - Ensures file exists on Windows, macOS, Linux.
    - Users can pass, or similar:
        '--image_icc_profile sRGB.icm' (current folder)
        '--image_icc_profile /System/Library/ColorSync/Profiles/sRGB Profile.icc'
        '--image_icc_profile ./profiles/AdobeRGB1998.icc'
        '--image_icc_profile ../Color/MyProfile.icc' (one folder level up)
        '--image_icc_profile C:/Color/MyProfile.icm'
        '--image_icc_profile /usr/share/color/icc/DisplayP3.icc'

    profile_arg : string provided to --pre_cond_profile
    arg_name    : the argument name for error messages

    Returns absolute validated path.
    Raises FileNotFoundError if missing or unreadable.
    """
    # --------------------------------------------------
    # CASE 1: Argument not provided → use fallback
    # --------------------------------------------------
    if profile_arg is None:
        return "sRGB.icm"

    # --------------------------------------------------
    # CASE 2: Argument provided but empty → ERROR
    # Example: --pre_cond_profile ""  or   --pre_cond_profile <space>
    # --------------------------------------------------
    if profile_arg.strip() == "":
        raise FileNotFoundError(
            f"Error: argument '{arg_name}' was provided but no file path was given."
        )

    profile_arg = profile_arg.strip()

    # Determine folder where script is located
    script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

    # --------------------------------------------------
    # CASE 3: Path-like input → check directly
    # --------------------------------------------------
    if os.path.sep in profile_arg or profile_arg.startswith("."):
        candidate = os.path.abspath(profile_arg)
        if os.path.isfile(candidate):
            return candidate
        else:
            raise FileNotFoundError(
                f"Error: ICC/ICM profile provided to '{arg_name}' not found:\n"
                f"  {profile_arg}\n"
                f"Please provide a valid existing file path."
            )

    # --------------------------------------------------
    # CASE 4: Bare filename → search script folder
    # --------------------------------------------------
    candidate = os.path.join(script_dir, profile_arg)
    if os.path.isfile(candidate):
        return os.path.abspath(candidate)

    raise FileNotFoundError(
        f"Error: ICC/ICM profile '{profile_arg}' for '{arg_name}' was not found.\n"
        f"Search attempted:\n"
        f"  - Provided path: {profile_arg}\n"
        f"  - Script folder: {script_dir}\n"
        f"Please provide a correct file path.\n"
        f"For simplicity the ICC/ICM profile can be copied to the folder location of this script."
    )


def resolve_output_paths(args, image_basename):
    """
    Resolve and validate --output argument.
    Returns (out_dir, out_base) where:
      - out_dir: absolute directory path (will be created later)
      - out_base: base filename without extension
    Raises SystemExit on error.
    """
    if not args.output:
        return os.getcwd(), image_basename

    out_path = os.path.abspath(os.path.expanduser(args.output))

    if os.path.isdir(out_path):
        return out_path, image_basename
    elif not os.path.exists(out_path) and args.output.endswith(os.sep):
        print(f"Error: Directory '{out_path}' does not exist.")
        sys.exit(1)
    else:
        # File path prefix
        out_dir = os.path.dirname(out_path) or os.getcwd()
        out_base = os.path.splitext(os.path.basename(out_path))[0]
        return out_dir, out_base


def Count_White(patches, white_point_xyz):
    """Count patches that are white (within tolerance of white point)."""
    if white_point_xyz is None:
        return 0

    tolerance = 0.005  # 0.005% tolerance for white detection
    white_count = 0

    for p in patches:
        if p.xyz100[1] is not None and abs(p.xyz100[1] - white_point_xyz[1]) <= tolerance:
            white_count += 1

    return white_count


def Count_Black(patches, black_point_xyz):
    """Count patches that are black (within tolerance of black point)."""
    if black_point_xyz is None:
        return 0

    tolerance = 0.005  # 0.005% tolerance for black detection
    black_count = 0

    for p in patches:
        if p.xyz100[1] is not None and abs(p.xyz100[1] - black_point_xyz[1]) <= tolerance:
            black_count += 1

    return black_count


def CountCompGrey(patches):
    """Count composite gray patches (equal device channel values R=G=B).

    Returns count of patches with equal RGB channel values.
    This corresponds to COMP_GREY_STEPS in ArgyllCMS targen specification.

    Args:
        patches: List of PatchInfo objects with rgb_percent attribute

    Returns:
        int: Number of composite gray patches
    """
    comp_grey_count = 0
    rgb_tolerance = 0.1   # Very tight RGB tolerance for balance
    min_lightness = 1.0   # Minimum RGB value for gray (above black)
    max_lightness = 99.0  # Maximum RGB value for gray (below white)

    for p in patches:
        # Check RGB balance (R=G=B within tolerance)
        if p.rgb_percent:
            r, g, b = p.rgb_percent

            # Check if all RGB channels are balanced
            if (abs(r - g) <= rgb_tolerance and
                abs(r - b) <= rgb_tolerance and
                abs(g - b) <= rgb_tolerance):
                # Also ensure it's not pure black or white
                if min_lightness <= r <= max_lightness:
                    comp_grey_count += 1

    return comp_grey_count


def CountSingleChannelRamps(patches, num_cols, num_rows):
    """Count single-dimension color ramps where two RGB channels remain static while the third creates a ramp pattern.

    Returns the most commonly detected number of steps for single-channel ramps.
    A single-channel ramp is defined as a sequence where:
    - Two RGB channels remain constant (within tolerance)
    - One RGB channel changes progressively (positive or negative direction)
    - At least 3 patches are required to form a valid ramp
    - Values are approximately evenly spaced (uniform spacing between steps)
    - Patches can be in any order (randomized)

    Args:
        patches: List of PatchInfo objects with rgb_percent attribute
        num_cols: Number of columns in the grid (unused, kept for compatibility)
        num_rows: Number of rows in the grid (unused, kept for compatibility)

    Returns:
        int: Most commonly detected number of steps for single-channel ramps (0 if no ramps found)
    """
    if not patches:
        return 0

    # Filter patches with valid RGB data
    rgb_patches = [p for p in patches if p.rgb_percent]
    if len(rgb_patches) < 3:
        return 0

    tolerance = 0.5  # Tolerance for detecting "static" channels
    min_ramp_length = 3  # Minimum patches to form a ramp
    ramp_step_counts = []  # Store step counts for all detected ramps

    # Extract RGB values
    rgb_values = [(p.rgb_percent[0], p.rgb_percent[1], p.rgb_percent[2]) for p in rgb_patches]

    def check_uniform_spacing(values):
        """Check if values are approximately evenly spaced."""
        if len(values) < 2:
            return True

        # Calculate spacing between consecutive values
        spacings = [values[i+1] - values[i] for i in range(len(values)-1)]

        # Check if all spacings are approximately equal
        if len(spacings) < 2:
            return True

        avg_spacing = sum(spacings) / len(spacings)
        spacing_tolerance = avg_spacing * 0.1  # 10% tolerance

        for spacing in spacings:
            if abs(spacing - avg_spacing) > spacing_tolerance:
                return False

        return True

    # Check each possible ramp pattern (R-static, G-static, B-static)
    for channel_idx, channel_name in [(0, 'R'), (1, 'G'), (2, 'B')]:
        # Get indices of the two static channels
        static_channels = [i for i in [0, 1, 2] if i != channel_idx]

        # Group patches by similar static channel values
        ramp_groups = {}

        for i, (r, g, b) in enumerate(rgb_values):
            # Get values of static channels
            static_key = (round(rgb_values[i][static_channels[0]] / tolerance) * tolerance,
                        round(rgb_values[i][static_channels[1]] / tolerance) * tolerance)

            if static_key not in ramp_groups:
                ramp_groups[static_key] = []
            ramp_groups[static_key].append((i, rgb_values[i][channel_idx]))

        # Analyze each group for uniform spacing ramp patterns
        for static_key, channel_values in ramp_groups.items():
            if len(channel_values) < min_ramp_length:
                continue

            # Sort by the changing channel value
            channel_values.sort(key=lambda x: x[1])

            # Extract just the channel values for analysis
            values_only = [val for idx, val in channel_values]

            # Check for monotonic increase or decrease
            is_increasing = all(values_only[i] <= values_only[i+1] for i in range(len(values_only)-1))
            is_decreasing = all(values_only[i] >= values_only[i+1] for i in range(len(values_only)-1))

            if is_increasing or is_decreasing:
                # Check if values are approximately evenly spaced
                if check_uniform_spacing(values_only):
                    ramp_steps = len(values_only)
                    ramp_step_counts.append(ramp_steps)

                    # DEBUG: Print ramp details when debug flag is enabled
                    if DEBUG:
                        static_channel_names = ['R', 'G', 'B']
                        print(f"\n=== DETECTED {channel_name}-CHANNEL RAMP ===")
                        print(f"Static channels: {static_channel_names[static_channels[0]]}={static_key[0]:.1f}%, {static_channel_names[static_channels[1]]}={static_key[1]:.1f}%")
                        print(f"Ramp direction: {'INCREASING' if is_increasing else 'DECREASING'}")
                        print(f"Steps: {ramp_steps}")
                        print(f"Channel values: {[f'{v:.1f}%' for v in values_only]}")
                        print("Patch details (SAMPLE_ID, RGB):")
                        for idx, val in channel_values:
                            patch = rgb_patches[idx]
                            r, g, b = patch.rgb_percent
                            print(f"  SAMPLE_ID {patch.index}: RGB=({r:.6f}, {g:.6f}, {b:.6f})")

    # Return the most commonly detected step count
    if not ramp_step_counts:
        return 0

    # Count frequency of each step count
    from collections import Counter
    step_frequencies = Counter(ramp_step_counts)

    # Return the most common step count (if tie, returns the first encountered)
    most_common_steps = step_frequencies.most_common(1)[0][0]
    return most_common_steps


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
        col_labels: List[str],
        density_extreme_records=None,
        device_combination_records=None
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
        self.density_extreme_records = density_extreme_records or []
        self.device_combination_records = device_combination_records or []

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
            fh.write('ORIGINATOR "read_image_patch_colors.py"\n')
            fh.write(f'CREATED "{created_ts}"\n')
            fh.write('ACCURATE_EXPECTED_VALUES "true"\n')

            # Write approximate white point if available
            if self.white_point_xyz is not None:
                Xw, Yw, Zw = [float(v) for v in self.white_point_xyz]
                fh.write(f'APPROX_WHITE_POINT "{Xw:.6f} {Yw:.6f} {Zw:.6f}"\n')
            else:
                fh.write('APPROX_WHITE_POINT "0.000000 0.000000 0.000000"\n')

            fh.write('COLOR_REP "iRGB"\n')

            # Add white, black, composite grey, and neutral patch counts
            white_count = Count_White(self.patches, self.white_point_xyz)
            black_count = Count_Black(self.patches, self.black_point_xyz)
            comp_grey_count = CountCompGrey(self.patches)
            single_channel_ramp_count = CountSingleChannelRamps(self.patches, self.num_cols, self.num_rows)

            # Only write count tags when values are non-zero
            if white_count > 0:
                fh.write(f'WHITE_COLOR_PATCHES "{white_count}"\n')
            if black_count > 0:
                fh.write(f'BLACK_COLOR_PATCHES "{black_count}"\n')
            if comp_grey_count > 0:
                fh.write(f'COMP_GREY_STEPS "{comp_grey_count}"\n')
            if single_channel_ramp_count > 0:
                fh.write(f'SINGLE_DIM_STEPS "{single_channel_ramp_count}"\n')
            fh.write('\n')

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
            fh.write('END_DATA\n\n')

            # =========================================
            # Dynamic DENSITY_EXTREME_VALUES block
            # =========================================
            fh.write("CTI1\n\n")
            fh.write('DESCRIPTOR "Argyll Calibration Target chart information 1"\n')
            fh.write('ORIGINATOR "Argyll targen"\n')
            fh.write(f'DENSITY_EXTREME_VALUES "{len(self.density_extreme_records)}"\n')
            fh.write(f'CREATED "{created_ts}"\n\n')
            # INDEX + RGB(3) + XYZ(3) = 7 fields
            fh.write('NUMBER_OF_FIELDS "7"\n')
            fh.write('BEGIN_DATA_FORMAT\n')
            fh.write("INDEX RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z\n")
            fh.write("END_DATA_FORMAT\n\n")
            fh.write(f'NUMBER_OF_SETS "{len(self.density_extreme_records)}"\n')
            fh.write('BEGIN_DATA\n')
            for idx, rgb, xyz, lab in self.density_extreme_records:
                R, G, B = rgb
                X, Y, Z = xyz
                fh.write(
                    f"{idx} "
                    f"{R:.6f} {G:.6f} {B:.6f} "
                    f"{X:.6f} {Y:.6f} {Z:.6f}\n"
                )
            fh.write('END_DATA\n\n')

            # =========================================
            # Dynamic DEVICE_COMBINATION_VALUES block
            # =========================================
            fh.write("CTI1\n\n")
            fh.write('DESCRIPTOR "Argyll Calibration Target chart information 1"\n')
            fh.write('ORIGINATOR "Argyll targen"\n')
            fh.write(f'DEVICE_COMBINATION_VALUES "{len(self.device_combination_records)}"\n')
            fh.write(f'CREATED "{created_ts}"\n\n')
            # INDEX + RGB(3) + XYZ(3) = 7 fields
            fh.write('NUMBER_OF_FIELDS "7"\n')
            fh.write('BEGIN_DATA_FORMAT\n')
            fh.write("INDEX RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z\n")
            fh.write("END_DATA_FORMAT\n\n")
            fh.write(f'NUMBER_OF_SETS "{len(self.device_combination_records)}"\n')
            fh.write('BEGIN_DATA\n')
            for idx, rgb, xyz, lab in self.device_combination_records:
                R, G, B = rgb
                X, Y, Z = xyz
                fh.write(
                    f"{idx} "
                    f"{R:.6f} {G:.6f} {B:.6f} "
                    f"{X:.6f} {Y:.6f} {Z:.6f}\n"
                )
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
            fh.write('ORIGINATOR "read_image_patch_colors.py"\n')
            fh.write(f'CREATED "{created_ts}"\n')
            fh.write('ACCURATE_EXPECTED_VALUES "true"\n')

            # Write approximate white point if available
            if self.white_point_xyz is not None:
                Xw, Yw, Zw = [float(v) for v in self.white_point_xyz]
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
    # --- Validate image path early ---
    if not args.image or not isinstance(args.image, str):
        print("Error: --image path not provided.")
        sys.exit(1)

    # Normalize (handles Windows/macOS/Linux)
    img_path = os.path.abspath(os.path.expanduser(args.image))

    # Basic extension check (optional but helpful)
    valid_exts = {".tif", ".tiff", ".jpg", ".jpeg", ".png", ".bmp"}
    ext = os.path.splitext(img_path)[1].lower()

    if ext not in valid_exts:
        print(f"Warning: Image file extension '{ext}' is uncommon or unsupported.")
        print("Supported extensions include:", ", ".join(sorted(valid_exts)))
        # Not fatal—continue.

    # Now check file existence
    if not os.path.isfile(img_path):
        print("Error: Image file does not exist:")
        print("  ", img_path)
        sys.exit(1)

    # Overwrite args.image with resolved absolute path
    args.image = img_path

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

    # Resolve and validate output paths early
    base_out = os.path.splitext(os.path.basename(args.image))[0]
    out_dir, out_base = resolve_output_paths(args, base_out)

    # Load image and compute geometry
    img16, _ = load_image_as_16bit_rgb(args.image)
    pre_cond_profile = args.pre_cond_profile
    #only do if pre_cond_profile exists
    if pre_cond_profile:
        pre_cond_profile = pre_cond_profile.lower()

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
        pre_cond_profile,
        args.sample_mode,
        args.output_order,
        args.mirror_output,
        args.rotate_grid
    )

    # Build data for Density Extreme Values in TI1 file
    density_extreme_records = build_patch_records(
        density_extr_val_list,
        pre_cond_profile
    )
    # Build data for Device Combination Values in TI1 file
    device_combination_records = build_patch_records(
        device_combi_val_list,
        pre_cond_profile
    )

    ti1_file = os.path.join(out_dir, f"{out_base}.ti1")
    ti2_file = os.path.join(out_dir, f"{out_base}.ti2")
    csv_file = os.path.join(out_dir, f"{out_base}.csv")

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
            density_extreme_records=density_extreme_records,
            device_combination_records=device_combination_records
        )
        writer.write()
        print("Wrote:", fname)

    print("Done.")


# ---------- CLI ----------
# Build argparse parser with expected command line options.
def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Read patch colours from grid image with label metadata.")
    p.add_argument('--image', '-i', required=True, help='Image file path')
    p.add_argument('--pre_cond_profile', default="sRGB.icm",
                   help='ICC/ICM profile file path, used as device pre-conditioning profile.')
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

    # --------------------------
    # Validate ICC/ICM profile
    # --------------------------
    try:
        args.pre_cond_profile = resolve_icc_profile_path(args.pre_cond_profile,"--pre_cond_profile")
    except FileNotFoundError as e:
        print(str(e))
        sys.exit(1)

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
