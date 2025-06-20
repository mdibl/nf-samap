#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-16
Version: 1.1.0
Purpose: Load SAM objects from a sample sheet CSV file and pickle them
"""

import argparse
import csv
from typing import NamedTuple
from pathlib import Path
from samalg import SAM
from log_utils import log
import pickle


class Args(NamedTuple):
    """Command-line arguments for the script"""

    sample_sheet: Path  # Path to the sample sheet CSV file
    output: Path # Path to the output directory


# --------------------------------------------------
def get_args() -> Args:
    """
    Parse command-line arguments.

    Returns:
        Args: NamedTuple containing parsed command-line arguments
    """

    parser = argparse.ArgumentParser(
        description="Load SAM objects from a sample sheet CSV file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--sample-sheet",
        metavar="FILE",
        type=Path,
        required=True,
        help="Path to the sample sheet CSV file containing 'id2' and 'h5ad' columns",
    )
    
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("."),
        help="Directory to save the SAM pickle outputs",
    )

    args = parser.parse_args()

    return Args(args.sample_sheet, args.output)


# --------------------------------------------------
def get_h5ad_dict(sample_sheet_path: Path) -> dict:
    """
    Read a sample sheet CSV file and return a dictionary mapping 'id2' to 'h5ad' file paths.

    Args:
        sample_sheet_path (Path): Path to the sample sheet CSV file.

    Returns:
        dict: A dictionary where the key is 'id2' and the value is the corresponding 'h5ad' file path.
    """
    h5ad_dict = {}
    with open(sample_sheet_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            h5ad_dict[row["id2"]] = row["h5ad"]
            log(f"  {row['id2']}: {row['h5ad']}", level="INFO")
    return h5ad_dict


# --------------------------------------------------
def load_sams(h5ad_dict: dict) -> dict:
    """
    Load SAM objects from a dictionary of h5ad file paths.

    Args:
        h5ad_dict (dict): A dictionary where the key is 'id2' and the value is the corresponding 'h5ad' file path.

    Returns:
        dict: A dictionary of SAM objects, keyed by 'id2'.
    """
    sams = {}
    for id2, h5ad in h5ad_dict.items():
        sams[id2] = SAM()
        sams[id2].load_data(h5ad)
        log(f"  Loading {id2}", level="INFO")
    return sams


# --------------------------------------------------
def pickle_sams(sams: dict, output_dir: Path) -> None:
    """
    Pickle each SAM object to a file named <id2>_sam.pkl in the specified output directory.

    Args:
        sams (dict): A dictionary of SAM objects, keyed by 'id2'.
        output_dir (Path): The directory where the pickled SAM objects will be saved.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for id2, sam in sams.items():
        out_path = output_dir / f"{id2}_sam.pkl"
        log(f"  Pickling {id2} to {out_path}", level="INFO")
        with open(out_path, "wb") as f:
            pickle.dump(sam, f)


# --------------------------------------------------
def main() -> None:
    """Load SAM objects from a sample sheet CSV file""" """
    Main function to load SAM objects from a sample sheet CSV file and pickle them.

    This function:
    1. Gets command-line arguments.
    2. Loads the sample sheet into a dictionary.
    3. Loads SAM objects from the dictionary.
    4. Pickles the SAM objects into the current directory.
    """
    
    log("Beginning execution of script", level="INFO")

    # Get command-line arguments
    args = get_args()
    log(f"Loaded sample sheet from {args.sample_sheet}", level="INFO")

    # Load the h5ad files as a dict from the sample sheet
    log("Loading h5ad dict", level="INFO")
    h5ad_dict = get_h5ad_dict(args.sample_sheet)
    log(f"Loaded h5ad dict with {len(h5ad_dict)} entries", level="INFO")

    # Load SAM objects from the h5ad dict
    log("Loading SAMs list", level="INFO")
    sams = load_sams(h5ad_dict)
    log(f"Loaded SAMs list with {len(sams)} entries", level="INFO")

    # Pickle SAM objects to files
    log("Pickling SAM objects", level="INFO")
    pickle_sams(sams, args.output)
    log(f"Script complete, see {args.output.resolve()}", level="INFO")


# --------------------------------------------------
if __name__ == "__main__":
    main()
