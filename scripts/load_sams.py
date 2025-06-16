#!/usr/bin/env python3
"""
Author : ryan <ryan@localhost>
Date   : 2025-06-16
Purpose: Load SAM objects from a sample sheet CSV file
"""

import argparse
import csv
from typing import NamedTuple
from pathlib import Path
from samalg import SAM
import pickle


class Args(NamedTuple):
    """ Command-line arguments """
    sample_sheet: Path

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Load SAM objects from a sample sheet CSV file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s',
                        '--sample-sheet',
                        metavar='FILE',
                        type=Path,
                        required=True,
                        help='Path to the sample sheet CSV file')

    args = parser.parse_args()

    return Args(args.sample_sheet)

# --------------------------------------------------
def get_h5ad_dict(sample_sheet_path: Path) -> dict:
    """Return a dict of {id2: h5ad} from the sample sheet CSV"""
    h5ad_dict = {}
    with open(sample_sheet_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            h5ad_dict[row['id2']] = row['h5ad']
    return h5ad_dict

# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    # Get command-line arguments
    args = get_args()
    print(f'Loaded sample_sheet: "{args.sample_sheet}"')
    
    # Load the h5ad files as a dict from the sample sheet
    h5ad_dict = get_h5ad_dict(args.sample_sheet)
    print(f'Loaded h5ad dict: {h5ad_dict}')
    
    # Load SAM objects from the h5ad dict
    sams = load_sams(h5ad_dict)
    print(f'Loaded {len(sams)} SAM objects')
    
    # Pickle SAM objects to files
    pickle_sams(sams, Path("."))

# --------------------------------------------------
def load_sams(h5ad_dict: dict) -> dict:
    """Load SAM objects from a dict of {id2: h5ad}"""
    sams = {}
    for id2, h5ad in h5ad_dict.items():
        sams[id2] = SAM()
        sams[id2].load_data(h5ad)
    return sams

# --------------------------------------------------
def pickle_sams(sams: dict, output_dir: Path) -> None:
    """Pickle each SAM object to a file named <id2>.pkl in the output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for id2, sam in sams.items():
        out_path = output_dir / f"{id2}.pkl"
        with open(out_path, "wb") as f:
            pickle.dump(sam, f)
        print(f"Pickled {id2} to {out_path}")

# --------------------------------------------------
if __name__ == '__main__':
    main()
