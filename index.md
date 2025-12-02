# read\_image\_patch\_colors
# Reference and Usage Guide
**Version:** 1.5<br>
**Author:** Knut Larsson<br>
**Purpose:** Generate ArgyllCMS compatible `.ti1` and `.ti2` files from a grid of color patches in an image.


## Table of Contents
- [Overview](#overview)
- [Use Cases](#use-cases)
- [Examples Included](#examples-included)
- [Command-line Arguments](#command-line-arguments)
- [Image Input](#image-input)
- [Grid Geometry](#grid-geometry)
- [Sampling](#sampling)
- [Labeling Rules](#labeling-rules)
- [Color Processing](#color-processing)
- [Output Files](#output-files)
- [Patch Value Ranges](#patch-value-ranges)
- [Valid Color Blocks](#valid-color-blocks)
- [Dependencies](#dependencies)
- [Notes on Accuracy](#notes-on-accuracy)
- [Patch Sampling Method](#patch-sampling-method)
- [Output File Format (TI2)](#output-file-format-ti2)
- [Error Conditions](#error-conditions)
- [Workflow Summary](#workflow-summary)
- [Examples](#examples)


# Overview
This program extracts color values from a rectangular grid of color patches in an image, computes colorimetric values (RGB percentages, XYZ, Lab and approximate white point), applies row/column labeling rules, and writes three output files:

- ArgyllCMS .ti1 file
- ArgyllCMS .ti2 file
- CSV file (space-separated)

It supports sampling in mean or median mode, configurable patch geometry, numeric or alphabetic labeling patterns, and RGB/XYZ/Lab output combinations.

The script is designed for use in color-management workflows where printed color targets are scanned or photographed and converted into Argyll measurement files. Alternatively, convert original calibration target images to attain `.ti1`, `.ti2` and `.csv` files with their patch color data.


# Use Cases
Example use cases for reading the image of a reference target:

1. Create a `.ti1` file so that one may use ArgyllCMS `printtarg` command and generate a target using the colors of the image, which then can be used with ArgyllCMS `chartread` and `colprof` to create a printer profile.<br>
	- Workflow:<br>
Read image with script → input `.ti1` to `printtarg` → Print generated Target image without color management → input `.ti2` to `chartread` then read target paches with colorimeter → input `.ti3` to `colprof` → output profile.

2. Create a `.ti2` file so that one may use ArgyllCMS `chartread`, and use a colorimeter to read the color values using the original SpyderPrint Targets, which then can be used with ArgyllCMS `colprof` to create a printer profile.
	- Workflow:<br>
Read image with script → Print original SpyderPrint Target image without color management → input `.ti2` to `chartread` then read target paches with colorimeter → input `.ti3` to `colprof` → output profile.

3. Create a `.ti1` file so that one may use ArgyllCMS `fakeread` command using colors of the image, which then can be used with ArgyllCMS `colprof` to create a simulated profile.<br>
	- Workflow:<br>
Read image with script  → input `.ti1` to `fakeread` → input `.ti3` to `colprof` → output simulated profile.


# Examples Included
See chapter [Examples](#examples) for detailed examples on how to read patches from different size patch grid images. All examples are based on Datacolor SpyderPrint Targets, which are:

- Expert Target (3-pages, 729-patches)
- Expert Target (large) (1-page, 729-patches)
- Expert Target Plus Grays (4-pages, 967-patches)
- EZ 729 Colors Plus Grays (9-pages, 996-patches)
- High Quality Target (1-page, 225-patches)
- High Quality Target Plus Grays (2-pages, 463-patches)

Target images used for creating .ti1, .ti2 and .csv files are under folder [Example Targets Read](https://github.com/soul-traveller/read_image_patch_colors/tree/main/Example%20Targets%20Read).

Several examples of generating an ArgyllCMS printtarg targets for SpyderPrint Targets are also included.


# Command-line Arguments
<table>
  <thead>
    <tr>
      <th>Argument</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>--image</code> / <code>-i</code></td>
      <td><b>(Required)</b> Path to the input image containing the colour-patch grid. Image must be a RGB image used as target for calibration.</td>
    </tr>

    <tr>
      <td><code>--pre_cond_profile</code></td>
      <td>
        ICC/ICM profile file path, used as device pre-conditioning profile. Default: 'sRGB.icm'. This adapts/shapes the XYZ and LAB values to a known device/printer profile or color space profile. This also affects the approx. white point value. For simplicity, place icc file in same folder as script.
        <br /><br />
        Supported path formats:
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;">• MyProfile.icm</td><td style="border: none;">(current folder)</td></tr>
          <tr><td style="border: none;">• /System/Library/ColorSync/Profiles/MyProfile.icc</td><td style="border: none;"></td></tr>
          <tr><td style="border: none;">• ./profiles/AdobeRGB1998.icc</td><td style="border: none;">(Current folder is ref.)</td></tr>
          <tr><td style="border: none;">• ../MyProfile.icc</td><td style="border: none;">(one folder up)</td></tr>
          <tr><td style="border: none;">• C:\Windows\MyProfile.icm</td><td style="border: none;"></td></tr>
          <tr><td style="border: none;">• C:/Color/MyProfile.icm</td><td style="border: none;"></td></tr>
          <tr><td style="border: none;">• /usr/share/color/icc/MyProfile.icc</td><td style="border: none;"></td></tr>
        </table>
      </td>
    </tr>

    <tr>
      <td><code>--patch_first_xy</code></td>
      <td><b>(Required)</b> "X,Y" coordinates of the <b>centre of the first patch</b> (top-left). Accepts integers or floats.</td>
    </tr>

    <tr>
      <td><code>--patch_last_xy</code></td>
      <td><b>(Required)</b> "X,Y" coordinates of the <b>centre of the last patch</b> (bottom-right). Accepts integers or floats.<br /><b>Note:</b> In some cases some pathces on last row of image may not exist. In these cases to ccordinates to the last row and column, as if it would exist, must be provided. After generating of .ti1, .ti2 and .csv files the last patches must be manually removed, if output must be exact like original image. If patch_label_order is col_then_row then the patces to remove are in sequence at the end of the list, but if row_then_col is used, then you would need to search for the patch lables (SAMPLE_LOC) to locate them and remove them.</td>
    </tr>

    <tr>
      <td><code>--patch_width_height_ratio</code></td>
      <td><b>(Required)</b> Width ÷ Height ratio (W/H) of a single patch.<br />Examples: <code>1.0</code> (square), <code>1.378</code>, etc.</td>
    </tr>

    <tr>
      <td><code>--num_cols</code></td>
      <td><b>(Required)</b> Number of columns in the grid.</td>
    </tr>

    <tr>
      <td><code>--num_rows</code></td>
      <td><b>(Required)</b> Number of rows in the grid.</td>
    </tr>

    <tr>
      <td><code>--sample_fraction</code></td>
      <td>
        Fraction of patch area to sample (float). Default: <b>0.50</b> (20%).<br />Constraints: <code>0 &lt; f ≤ 0.6</code>.<br />The sampling area is a centered square clamped to at least <b>3×3 px</b> and at most <b>60%</b> of patch size.
      </td>
    </tr>

    <tr>
      <td><code>--row_labels</code></td>
      <td>
        <b>(Required)</b> Label sequence for rows. Must match <code>num_rows</code> exactly.<br />
        Allowed patterns:
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;"><code>1-15</code></td><td style="border: none;">Numeric range</td></tr>
          <tr><td style="border: none;"><code>03-12</code></td><td style="border: none;">Numeric range (zero-padding preserved)</td></tr>
          <tr><td style="border: none;"><code>0001-0120</code></td><td style="border: none;">Numeric range (zero-padding preserved)</td></tr>
          <tr><td style="border: none;"><code>A-Z</code></td><td style="border: none;">Alphabetic range</td></tr>
          <tr><td style="border: none;"><code>A-AC</code></td><td style="border: none;">Alphabetic range</td></tr>
          <tr><td style="border: none;"><code>BQ-CF</code></td><td style="border: none;">Alphabetic range</td></tr>
        </table>
      </td>
    </tr>

    <tr>
      <td><code>--col_labels</code></td>
      <td><b>(Required)</b> Label sequence for columns (same rules as <code>row_labels</code>).</td>
    </tr>

    <tr>
      <td><code>--patch_label_order</code></td>
      <td>
        <b>(Required)</b> Determines patch label composition:<br />
        (Examples assume row_labels are numeric and col_labels are alphabetic)
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;"><code>col_then_row</code></td><td style="border: none;">→ e.g. <code>A12</code></td></tr>
          <tr><td style="border: none;"><code>row_then_col</code></td><td style="border: none;">→ e.g. <code>12A</code></td></tr>
        </table>
      </td>
    </tr>

    <tr>
      <td><code>--output_color_space</code></td>
      <td>
        <b>(Required)</b> Comma-separated list specifying which colour spaces to output (in order). Allowed tokens:<br />
        • <code>RGB</code><br />
        • <code>XYZ</code><br />
        • <code>LAB</code><br />
        Rules:<br />
        • TI1 can include only RGB and/or XYZ<br />
        • TI2/CSV include all listed tokens<br />
        Examples:
			<table style="border-collapse: collapse; border: none;">
			  <tr>
			    <td style="border: none;"><code>RGB,XYZ</code></td>
			    <td style="border: none;">→ TI1: RGB,XYZ; TI2: RGB,XYZ.</td>
			  </tr>
			  <tr>
			    <td style="border: none;"><code>XYZ,LAB,RGB</code></td>
			    <td style="border: none;">→ TI1: XYZ,RGB; TI2: XYZ,LAB,RGB.</td>
			  </tr>
			</table>
      </td>
    </tr>

    <tr>
      <td><code>--sample_mode</code></td>
      <td>
        Patch sampling aggregation method:<br />
			<table style="border-collapse: collapse; border: none;">
			  <tr>
			    <td style="border: none;"><code>mad <i>(default)</i></code></td>
			    <td style="border: none;">→ robust mean via MAD-based sigma clipping</td>
			  </tr>
			  <tr>
			    <td style="border: none;"><code>mean</code></td>
			    <td style="border: none;">→ arithmetic mean (not robust)</td>
			  </tr>
			  <tr>
			    <td style="border: none;"><code>median</code></td>
			    <td style="border: none;">→ robust, but biases asymmetric patches</td>
			  </tr>
			</table>
      </td>
    </tr>

    <tr>
      <td><code>--output_order</code></td>
      <td>
        Determines initial traversal order before rotation/mirroring, if applied:
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;"><code>row_major</code></td><td style="border: none;">iterate row by row (top→bottom, left→right)</td></tr>
          <tr><td style="border: none;"><code>column_major</code></td><td style="border: none;">iterate column by column (left→right, top→bottom)</td></tr>
        </table>
      </td>
    </tr>

    <tr>
      <td><code>--mirror_output</code></td>
      <td>
        Reverse traversal <b>after</b> <code>output_order</code> and <code>rotate_grid</code>:
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;"><code>row_major</code></td><td style="border: none;">reverse column labels for each row (vertical flip)</td></tr>
          <tr><td style="border: none;"><code>column_major</code></td><td style="border: none;">reverse row labels for each column (horizontal flip)</td></tr>
        </table>
      </td>
    </tr>

    <tr>
      <td><code>--rotate_grid</code></td>
      <td>
        Rotate the entire patch grid <b>clockwise</b> before mirroring (if <code>mirror_output</code> applied). Allowed values:
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;"><code>0</code></td><td style="border: none;">(default)</td></tr>
          <tr><td style="border: none;"><code>90</code></td><td style="border: none;">rotate 90° CW</td></tr>
          <tr><td style="border: none;"><code>180</code></td><td style="border: none;">rotate 180°</td></tr>
          <tr><td style="border: none;"><code>270</code></td><td style="border: none;">rotate 270° CW</td></tr>
        </table>
        Applied <b>after output_order</b>, <b>before mirror_output</b>.
      </td>
    </tr>

    <tr>
      <td><code>--output</code></td>
      <td>
        Base filename and/or path for generated <code>.ti1</code>, <code>.ti2</code>, and <code>.csv</code> output files.<br />
        If omitted, filenames are based on the input image name and placed in the same folder as the script. If path is specified, the script will place all output files. If the specified path also includes a filename, the script will use that filename for all output files.
        <table style="border-collapse: collapse; border: none;">
          <tr><td style="border: none;">--output "./Outputfolder/"</td><td style="border: none;">Only directory, current folder is reference (keep input image name)</td></tr>
          <tr><td style="border: none;">--output "./Outputfolder/Outputfilename"</td><td style="border: none;">Directory + filename, current folder is reference (override input image name)</td></tr>
          <tr><td style="border: none;">--output "User/Username/ Outputfolder/Outputfilename"</td><td style="border: none;">Directory + filename (override input image name)</td></tr>
          <tr><td style="border: none;">--output "../Outputfolder/"</td><td style="border: none;">Directory in parallel folder, one up in structure (keep input image name)</td></tr>
        </table>
      </td>
    </tr>

    <tr>
      <td><code>--debug</code></td>
      <td>Enable diagnostic printing for the first few patches. Shows sampling geometry, pixel statistics, and intermediate RGB calculations.</td>
    </tr>

  </tbody>
</table>



# Image Input
Accepted image types: any PIL-compatible raster file
Color depth: handled as 16-bit RGB internally

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

Sampling modes (parameter `--sample_mode`):

- mad    – robust mean using MAD-based sigma clipping (default)
- mean   – plain arithmetic mean (sensitive to outliers)
- median – pure median (very robust, but may bias bright/dark values)

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

Color conversions (XYZ, Lab) are computed through ArgyllCMS xicclu command with 
selected color icc profile (default `sRGB.icm`). Absolute colorimetric intent with D50 illuminant is used, same as how ArgyllCMS targen creates patch colors, and is expected by printtarg and chartread. Targen uses a sRGB-like model for generating XYZ and LAB D50, unless another profile i provided. This is why `sRGB.icm` color space profile is used by default, but may be overridden by choosing a different profile with `--pre_cond_profile`.

XYZ output:

- CIEXYZ (CIE 1931 2° standard observer)
- D50 white reference (Profile Connection Space in ICC)
- Absolute colorimetric intent
- White scaled to Y=100, unbounded positive, values 0-100
- Numeric precision: six decimal places.

Lab output:

- CIE L* a* b* (1976), computed using **D50** reference white (Xn, Yn, Zn)
- Absolute colorimetric intent
- L in range ~0..100, a and b range typical -128..128.
- Numeric precision: six decimal places.

- All colorimetric output is **D50** and **linear** (no gamma applied in XYZ or Lab).

The patch with the highest Y (XYZ_Y) is stored as APPROX\_WHITE\_POINT. 
Note that XYZ Absolute colorimetric D50 is almost identical to XYZ Relative colorimetric D65, which can confuse some users to think that the APPROX\_WHITE\_POINT is computed using D65. This is not the case.


# Output Files
Three files:

### 1. TI1 (.ti1)
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

#### COMP\_GREY\_STEPS
Number represents count of patches with equal RGB channel values, basically neutral colors. Note: This is not the NEUTRAL_STEPS tag, as they would be generated by ArgyllCMS targen taking into acount the neutral axis for a given icc profile. This script does not support creating NEUTRAL_STEPS count.

RGB Balance Check:<br>

- RGB Tolerance: 0.1%. Very tight tolerance requiring near-perfect R=G=B balance
- Lightness Range: 1.0% to 99.0%. Excludes pure black and white patches.
- Logic: A patch is neutral if all RGB channels are balanced within 0.1% and 
  fall within the specified lightness range.

#### SINGLE\_DIM\_STEPS
This number is estimated by reading single-dimension color ramps where two RGB channels remain static while the third creates a ramp pattern. The number represents the most commonly detected number of steps for single-channel ramps. A single-channel ramp is defined as a sequence where:<br>

  - Two RGB channels remain constant (within tolerance)
  - One RGB channel changes progressively (positive or negative direction)
  - At least 3 patches are required to form a valid ramp
  - Values are approximately evenly spaced (uniform spacing between steps)
  - Patches can be in any order (randomized)

<b>Note: if merging two or more `.ti1` files the value of this field is not added like other tags. Instead, use the largest number in the merged file.</b>

### 2. TI2 (.ti2)
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

SAMPLE\_LOC corresponds to labels such as "A1", "D14", "C07", etc.

### 3. CSV (.csv)
Space-separated table with SAMPLE\_ID, SAMPLE\_LOC, and selected color values.


# Patch Value Ranges
RGB: 0–100<br>
XYZ: typically 0–100<br>
Lab_L: 0–100<br>
Lab_a/b: approx. –128 to +128<br>

All values have six decimal places.


# Valid Color Blocks
`--output_color_space` must include at least one of: RGB, XYZ, LAB

Example:
`--output_color_space RGB,XYZ`

Order controls column order in TI1/TI2/CSV.


# Dependencies
- Python 3.8+
- numpy
- Pillow
- ArgyllCMS xicclu


# Notes on Accuracy
- This script seeks to minimize rounding/quantization error:
  * If the image has 16-bit/channel data it will be used directly.
  * If the image is 8-bit/channel (most typical images), it will be scaled up
    exactly to 16-bit (value * 257) before colour conversion to reduce rounding
    / mapping error in conversions.
  * Conversion from RGB -> XYZ -> Lab is done using ArgyllCMS xicclu with
  precision floats; XYZ is scaled so Yn = 100 and written with 6 decimal places.

- This script produces **colorimetrically correct** PCS values
  for the image that was measured.


# Patch Sampling Method
For each patch:

1. Patch centre (Cx, Cy) is computed by linear interpolation between
   patch_first_xy → patch_last_xy across the grid.
2. A square sampling region of size `sample_size × sample_size`
   (odd number, min 3, max 60% of patch) is centred at (Cx, Cy).
3. The script extracts all pixels inside this region.
4. The sampled RGB values are averaged (MAD-based sigma) in 16-bit integer 
	space (8-bit inputs are scaled to 16-bit via *257).
5. Colour conversions are applied using xicclu:
      - RGB -> CIEXYZ XYZ (D50) 
      - RGB -> Lab (D50)


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
Use `--debug` flag to investigate issues relating to color on patches.

Script stops if:

- Coordinates and Row/Column numbers cause measurement areas to be missaligned with patch centres of the target image.
- Labels do not match row/column counts
- Output color space tokens are invalid
- sample_fraction exceeds allowed range
- Grid geometry cannot be computed
- Any I/O error occurs during file writing
- File path to image, pre_cond_profile, or output directory is invalid


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
  --patch_first_xy 72,51 \
  --patch_last_xy 1250,934 \
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
  --image "LaserSoft Advanced ISO12641-2 Target (864-patches).tif" \
  --patch_first_xy 168,168 \
  --patch_last_xy 1335,935 \
  --patch_width_height_ratio 1.0 \
  --num_cols 36 \
  --num_rows 24 \
  --row_labels A-X \
  --col_labels 1-36 \
  --patch_label_order row_then_col \
  --output_color_space RGB,XYZ
```
