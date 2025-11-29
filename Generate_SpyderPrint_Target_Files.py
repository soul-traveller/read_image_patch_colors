#!/usr/bin/env python3
"""
Generate_SpyderPrint_Target_Files.py
Version: 1.0
Author: Knut Larsson

================================================================================
OVERVIEW
--------
This script automates the generation of ArgyllCMS compatible target files (.ti1, .ti2, .csv)
from SpyderPrint target images. It runs read_image_patch_colors.py for each configured
target with the appropriate parameters, eliminating the need to manually type commands
for each image.

The script handles:
- Multiple target types and configurations
- Cross-platform execution (Windows, macOS, Linux)
- Proper path handling with spaces
- Consistent Python interpreter usage
- Batch processing of all configured targets

================================================================================
HOW TO USE
-----------
1. **Prerequisites:**
   - Ensure read_image_patch_colors.py is in the same directory
   - Install required packages: pip install numpy Pillow
   - Install ArgyllCMS (xicclu must be in PATH)
   - Place target images in "Example Targets Read/SpyderPrint Targets/"
   - Place ICC profiles in "Profiles/" (optional, includes sRGB.icm by default)

2. **Run the script:**
   ```bash
   python3 Generate_SpyderPrint_Target_Files.py
   ```

3. **Output:**
   - .ti1, .ti2, and .csv files are generated in each target's folder
   - Progress is displayed for each target being processed
   - Errors are reported if any target fails to process

================================================================================
ADDING NEW TARGETS
------------------
To add new target configurations, modify the `configurations` list in the [main()](cci:1://file:///Users/knut/Documents/Calibration/Scripts/ArgyllCMS%20scripts/Published%20on%20Github/read_image_patch_colors/read_image_patch_colors/Generate_SpyderPrint_Target_Files.py:11:0-40:60) function:

```
python
configurations = [
    # Existing configurations...

    # Add new target like this:
    ("Target Name", "ImageFilename.tif",
     "first_x,first_y", "last_x,last_y",
     "ratio", "cols", "rows",
     "row_labels", "col_labels"),
]
```

**Configuration Parameters:**
- Target Name: Descriptive name for the target (used in logging)
- Image Filename: Name of the image file in the target folder
- first_x,first_y: Center coordinates of top-left patch
- last_x,last_y: Center coordinates of bottom-right patch
- ratio: Width/Height ratio of individual patches
- cols: Number of columns in the grid
- rows: Number of rows in the grid
- row_labels: Label range for rows (e.g., "1-15", "A-D")
- col_labels: Label range for columns (e.g., "A-O", "1-12")

**Steps to add a new target:**
1. Place the target image in the appropriate folder
2. Measure the patch grid parameters (first/last coordinates, ratio, dimensions)
3. Determine the labeling scheme
4. Add the configuration tuple to the configurations list
5. Run the script to test the new configuration

**Example: Adding a new 12x8 target:**
```
python
("My Custom Target", "custom_12x8.tif",
 "100,80", "500,400", "1.2", "12", "8", "1-8", "A-L"),
```

================================================================================
FOLDER STRUCTURE
----------------
The script expects this folder structure:

```
read_image_patch_colors/
├── Generate_SpyderPrint_Target_Files.py
├── read_image_patch_colors.py
├── Example Targets Read/
│   └── SpyderPrint Targets/
│       ├── Expert Target Page 1 of 3 (243-patches).tif
│       ├── Expert Target (large)(729-patches).tif
│       └── ... (other target images)
├── Profiles/
│   ├── sRGB.icm
│   └── ... (other ICC profiles)
└── Output/ (generated files appear in target folders)
```

================================================================================
TROUBLESHOOTING
---------------
- **ModuleNotFoundError**: Ensure Pillow is installed and the same Python interpreter is used
- **xicclu not found**: Install ArgyllCMS and ensure it's in your PATH
- **File not found**: Check that target images exist in the expected folder
- **Coordinate errors**: Verify patch_first_xy and patch_last_xy match the actual grid
- **Label mismatches**: Ensure row_labels and col_labels match the grid dimensions

Use the --debug flag in read_image_patch_colors.py to investigate patch sampling issues.

================================================================================
CUSTOMIZATION
-------------
- Modify the default output color space by changing the --output_color_space parameter
- Adjust sampling mode by changing --sample_mode (mad, mean, median)
- Change sample_fraction to modify the sampling area within each patch
- Add additional ICC profiles to the Profiles folder and reference them in configurations

================================================================================
"""


import os
import sys
import subprocess
from pathlib import Path

def main():
    # Detect script directory and set paths
    script_dir = Path(__file__).parent
    targets_root = script_dir / "Example Targets Read/SpyderPrint Targets"
    profiles_root = script_dir / "Profiles"

    print("Processing SpyderPrint targets...")
    print(f"Target folder: {targets_root}")
    print(f"Profiles folder: {profiles_root}")
    print()

    # Configuration: (folder, image, first_xy, last_xy, ratio, cols, rows, row_labels, col_labels)
    configurations = [
        ("Expert Target (3-pages, 729-patches)", "Expert Target Page 1 of 3 (243-patches).tif",
         "48,45", "872,604", "1.075", "18", "14", "1-14", "A-R"),
        ("Expert Target (3-pages, 729-patches)", "Expert Target Page 2 of 3 (243-patches).tif",
         "48,45", "872,604", "1.075", "18", "14", "1-14", "A-R"),
        ("Expert Target (3-pages, 729-patches)", "Expert Target Page 3 of 3 (243-patches).tif",
         "48,45", "872,604", "1.075", "18", "14", "1-14", "A-R"),

        ("Expert Target (large) (1-page, 729-patches)", "Expert Target (large)(729-patches).tif",
         "72,51", "1250,934", "1.4", "27", "27", "1-27", "A-AA"),

        ("Expert Target Plus Grays (4-pages, 967-patches)", "Expert Target Plus Grays 1 of 4 (243-patches).tif",
         "48,45", "872,604", "1.075", "18", "14", "1-14", "A-R"),
        ("Expert Target Plus Grays (4-pages, 967-patches)", "Expert Target Plus Grays 2 of 4 (243-patches).tif",
         "48,45", "872,604", "1.075", "18", "14", "1-14", "A-R"),
        ("Expert Target Plus Grays (4-pages, 967-patches)", "Expert Target Plus Grays 3 of 4 (243-patches).tif",
         "48,45", "872,604", "1.075", "18", "14", "1-14", "A-R"),
        ("Expert Target Plus Grays (4-pages, 967-patches)", "Expert Target Plus Grays 4 of 4 (238-patches).tif",
         "56,51", "847,592", "1.09803921568", "17", "14", "1-14", "A-Q"),

        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 1 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 2 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 3 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 4 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 5 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 6 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 7 of 9 (108-patches).tif",
         "111,64", "1039,746", "1.89286", "9", "12", "1-12", "A-I"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 8 of 9 (120-patches).tif",
         "87,67", "1023,749", "1.7142857", "10", "12", "1-12", "A-J"),
        ("EZ 729 Colors Plus Grays (9-pages, 996-patches)", "EZ 729 Colors Plus Grays 9 of 9 (120-patches).tif",
         "87,67", "1023,749", "1.7142857", "10", "12", "1-12", "A-J"),

        ("High Quality Target (1-page, 225-patches)", "High Quality Target (225-patches).tif",
         "61,50", "848,612", "1.37837838", "15", "15", "1-15", "A-O"),

        ("High Quality Target Plus Grays (2-pages, 463-patches)", "High Quality Target Plus Grays 1 of 2 (225-patches).tif",
         "61,50", "848,612", "1.37837838", "15", "15", "1-15", "A-O"),
        ("High Quality Target Plus Grays (2-pages, 463-patches)", "High Quality Target Plus Grays 2 of 2 (238-patches).tif",
         "56,51", "847,592", "1.09803921568", "17", "14", "1-14", "A-Q"),
    ]

    # Process each configuration
    for config in configurations:
        process_config(targets_root, profiles_root, *config)

def process_config(targets_root, profiles_root, folder_name, image_name, first_xy, last_xy,
                  ratio, cols, rows, row_labels, col_labels):
    """Process a single configuration."""

    folder_path = targets_root / folder_name
    image_path = folder_path / image_name

    # Check if image exists
    if not image_path.exists():
        print(f"WARNING: Image not found: {image_path}")
        return

    print(f"Processing: {folder_name}/{image_name}")

    # Build command as a string with proper quoting
    cmd_parts = [
        sys.executable,  # Use the same Python interpreter
        "read_image_patch_colors.py",
        f"--image \"{image_path}\"",
        f"--pre_cond_profile \"{profiles_root / 'sRGB.icm'}\"",
        f"--patch_first_xy {first_xy}",
        f"--patch_last_xy {last_xy}",
        f"--patch_width_height_ratio {ratio}",
        f"--num_cols {cols}",
        f"--num_rows {rows}",
        f"--row_labels {row_labels}",
        f"--col_labels {col_labels}",
        "--patch_label_order col_then_row",
        "--output_color_space RGB,XYZ",
        f"--output \"{folder_path}\""
    ]

    # Create command string for shell execution
    cmd_string = " ".join(cmd_parts)

    # Print command
    print(f"    Command: {cmd_string}")

    # Execute the command with shell=True for proper quoting
    try:
        result = subprocess.run(cmd_string, shell=True, check=True, capture_output=True, text=True)
        print(f"    Success: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"    Error: {e}")
        if e.stderr:
            print(f"    Stderr: {e.stderr}")

    print()

if __name__ == "__main__":
    main()
