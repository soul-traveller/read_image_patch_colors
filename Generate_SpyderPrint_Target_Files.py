#!/usr/bin/env python3
"""
Cross-platform script to process SpyderPrint target images.
Generates and executes read_image_patch_colors.py commands for each target configuration.
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
