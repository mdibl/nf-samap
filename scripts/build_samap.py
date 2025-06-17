#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-16
Version: 1.0.0
Purpose: Build a SAMAP object from sams and f_maps
"""

import argparse
import csv
import pickle
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple
from pathlib import Path


class Args(NamedTuple):
    """ Command-line arguments for the script"""
    
    sams_dir: Path      # Directory containing SAM pickles
    sample_sheet: Path  # Path to the sample sheet CSV file
    f_maps: Path        # Path to the f_maps directory


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
        '-s', '--sample-sheet',
        required=True,
        type=Path,
        help='Path to the sample sheet CSV file'
    )

    parser.add_argument(
        '-m', '--f-maps',
        required=True,
        type=Path,
        help='Path to the f_maps directory'
    )

    args = parser.parse_args()
    return Args(args.sams_dir, args.sample_sheet, args.f_maps)


# --------------------------------------------------
def load_species_dict(sample_sheet_path: Path, sams_dir: Path) -> dict:
    """
    Load a dictionary of species, mapping id2 to corresponding SAM objects from the sams_dir directory.

    Args:
        sample_sheet_path (Path): Path to the sample sheet CSV file.
        sams_dir (Path): Path to the directory containing the SAM pickle files.

    Returns:
        dict: A dictionary with id2 as the key and the corresponding SAM object as the value.
    """
    species = {}
    with open(sample_sheet_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            id2 = row["id2"]
            # Find the pickle file in sams_dir that starts with id2 (first two characters)
            matching_files = list(sams_dir.glob(f"{id2}*.pkl"))
            if not matching_files:
                print(f"[WARN] No SAM pickle found for id2={id2} in {sams_dir}")
                continue
            sam_path = matching_files[0]
            with open(sam_path, "rb") as f:
                species[id2] = pickle.load(f)
            print(f"[INFO] Loaded SAM for {id2} from {sam_path}")
    return species


# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the species dictionary from the sample sheet and SAM files.
    3. Validates the f_maps directory.
    4. Creates a SAMAP object using the loaded species and f_maps data.
    5. Saves the SAMAP object to a pickle file.
    """

    # Parse command-line arguments
    args = get_args()
    print(f'[INFO] SAMs directory: {args.sams_dir}')
    print(f'[INFO] f_maps directory: {args.f_maps}')
    print(f'[INFO] Sample sheet: {args.sample_sheet}')
    
    # Load species dictionary from sample sheet
    species_dict = load_species_dict(args.sample_sheet, args.sams_dir)
    print(f'[INFO] Loaded species dictionary with {len(species_dict)} entries')
    for id2, sam_object in species_dict.items():
        print(f'  [INFO] {id2}: {sam_object}')
        
    # Check if f_maps directory exists and ensure compatibility with SAMap
    f_maps = str(args.f_maps)
    if not f_maps.endswith('/'):
        f_maps += '/'
    if not Path(f_maps).exists():
        raise FileNotFoundError(f"[ERROR] f_maps directory {f_maps} does not exist")
    else:
        print(f'[INFO] f_maps directory found: {f_maps}')
        for map_file in Path(f_maps).rglob('*.txt'):  # Use rglob for recursive search
            print(f'  [INFO] Found map file: {map_file}')

    # Create SAMAP object
    print('[INFO] Creating SAMAP object...')
    samap = SAMAP(
        sams=species_dict,
        f_maps=f_maps,
        save_processed=False,
    )
    print(f'[INFO] Created SAMAP object with {len(samap.sams)} SAMs')
    
    # Save SAMAP object
    print('[INFO] Saving SAMAP object...')
    save_samap(samap, str('samap.pkl'))
    print('[INFO] SAMAP object saved to samap.pkl')

# --------------------------------------------------
if __name__ == '__main__':
    main()
