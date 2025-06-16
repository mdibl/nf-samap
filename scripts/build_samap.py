#!/usr/bin/env python3
"""
Author : ryan <ryan@localhost>
Date   : 2025-06-16
Purpose: Build a SAMAP object from sams
"""

import argparse
import csv
import pickle
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple
from pathlib import Path


class Args(NamedTuple):
    """ Command-line arguments """
    sams_dir: Path
    sample_sheet: Path
    f_maps: Path


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

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
    """Load species dictionary from sample sheet CSV, mapping id2 to SAM objects from sams_dir"""
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
    """ Main entry point """

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
