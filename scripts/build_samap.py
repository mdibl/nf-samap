#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-16
Version: 1.0.0
Purpose: Build a SAMAP object from sams and maps
"""

import argparse
import csv
import pickle
import os
from log_utils import log
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple
from pathlib import Path


class Args(NamedTuple):
    """ Command-line arguments for the script"""
    
    sams_dir: Path      # Directory containing SAM pickles
    sample_sheet: Path  # Path to the sample sheet CSV file
    maps: Path          # Path to the maps directory
    name: str           # Name of the output pickle
    output_dir: Path    # Path to the output directory


# --------------------------------------------------
def get_args() -> Args:
    """
    Parse and return command-line arguments.

    Returns:
        Args: A named tuple containing parsed command-line arguments for sams_dir, sample_sheet, and f_maps.
    """
    parser = argparse.ArgumentParser(
        description='Build a SAMAP object from a directory of SAMs and a sample sheet',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-d', '--sams-dir',
        required=True,
        type=Path,
        help='Directory containing SAM pickle files'
    )

    parser.add_argument(
        '-i', '--id2',
        required=True,
        type=str,
        help='list of id2 from the Sample Sheet'
    )

    parser.add_argument(
        '-m', '--maps',
        required=True,
        type=Path,
        help='Path to the maps directory'
    )
    
    parser.add_argument(
        '-n', '--name',
        required=False,
        type=str,
        help='Name of the output pickle',
        default='samap.pkl'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    args = parser.parse_args()
    return Args(args.sams_dir, args.id2, args.maps, args.name, args.output_dir)


# --------------------------------------------------
def load_species_dict(id2: str, sams_dir: Path) -> dict:
    """
    Load a dictionary of species, mapping id2 to corresponding SAM objects from the sams_dir directory.

    Args:
        sample_sheet_path (Path): Path to the sample sheet CSV file.
        sams_dir (Path): Path to the directory containing the SAM pickle files.

    Returns:
        dict: A dictionary with id2 as the key and the corresponding SAM object as the value.
    """
    species = {}
    for val in id2
        log(f"  Attempting to load SAM pickle for '{val}'", "INFO")
        # Find the pickle file in sams_dir that starts with id2 
        matching_files = list(sams_dir.glob(f"{val}*.pkl"))
        if not matching_files:
            log(f"  No SAM pickle found for '{val}' in '{sams_dir}'", "ERROR")
        sam_path = matching_files[0]
        with open(sam_path, "rb") as f:
            species[val] = pickle.load(f)
        log(f"  Loaded SAM for '{val}' from '{sam_path}'", "INFO")
    return species


"""     with open(sample_sheet_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            id2 = row["id2"]
            log(f"  Attempting to load SAM pickle for '{id2}'", "INFO")
            # Find the pickle file in sams_dir that starts with id2 
            matching_files = list(sams_dir.glob(f"{id2}*.pkl"))
            if not matching_files:
                log(f"  No SAM pickle found for '{id2}' in '{sams_dir}'", "ERROR")
            sam_path = matching_files[0]
            with open(sam_path, "rb") as f:
                species[id2] = pickle.load(f)
            log(f"  Loaded SAM for '{id2}' from '{sam_path}'", "INFO")
    return species """


# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the species dictionary from the sample sheet and SAM files.
    3. Validates the maps directory.
    4. Creates a SAMAP object using the loaded species and maps data.
    5. Saves the SAMAP object to a pickle file.
    """

    # Parse command-line arguments
    log("Loading arguments", "INFO")
    args = get_args()
    sams_dir = args.sams_dir
    log(f"  Using SAMs directory '{sams_dir}'", "DEBUG")
    maps = str(args.maps)
    log(f"  Using maps directory '{maps}'", "DEBUG")
    id2 = args.id2
    log(f"  Using sample sheet '{id2}'", "DEBUG")
    name = args.name
    log(f"  SAMAP object will be saved with name '{name}'", "DEBUG")
    output_dir = args.output_dir
    log(f"  SAMAP object will be saved to '{output_dir}'", "DEBUG")
    
    # Load species dictionary from sample sheet
    log("Loading species dictionary from sample sheet", "INFO")
    species_dict = load_species_dict(id2, sams_dir)
    log(f"Loaded species dictionary with {len(species_dict)} entries", "INFO")


    # Ensure maps is valid and formatted correctly
    log(f"Ensuring validity of '{maps}'", "INFO")
    if not maps.endswith('/'): # SAMap *will* crash if passed a dir without a '/'
        maps += '/'
        log(f"Provided maps directory does not end with '/', changing to '{maps}'", "WARN")
    if not Path(maps).exists():
        error_message = f"Maps directory '{maps}' does not exist"
        log(error_message, "ERROR")
        raise FileNotFoundError(error_message)
    else:
        log(f"Maps directory found at '{maps}'", "INFO")
        for map_file in Path(maps).rglob('*.txt'):  # Use rglob for recursive search
            log(f"  Found map file '{map_file}", "DEBUG")


    # Create SAMAP object
    log("Attempting to create SAMAP object", "INFO")
    samap = SAMAP(
        sams=species_dict,
        f_maps=maps,
        save_processed=False,
    )
    log("Successfully created SAMAP object with {len(samap.sams)} SAMs", "INFO")
    
    # Save SAMAP object
    log("Attempting to pickle SAMAP object", "INFO")
    save_samap(samap, os.path.join(output_dir, name))
    log(f"Successfully pickled SAMAP object '{name}' to '{output_dir}'")

# --------------------------------------------------
if __name__ == '__main__':
    main()
