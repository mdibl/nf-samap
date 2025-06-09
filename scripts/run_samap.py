#!/root/miniconda/bin/python
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP with command-line input
"""

import argparse
import json
from typing import Dict
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple


class Args(NamedTuple):
    """ Command-line arguments """
    input: Path
    data: Path


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Run SAMAP with command-line input',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input',
                        type=Path,
                        default=Path('inputs.json'),
                        help='Path to JSON with species names and h5ad paths')

    parser.add_argument('--data',
                        type=Path,
                        default=Path('data/'),
                        help='Path to data directory containing maps and JSON file')

    args = parser.parse_args()

    return Args(args.input,
                args.data)


# --------------------------------------------------
def load_json(args: Args) -> Dict[str, Path]:
    """ Verify and load the JSON file with species names and h5ad paths """
    print(f"Loading JSON file: {args.input}")
    if not args.input.exists():
        raise ValueError(f"JSON file '{args.input}' does not exist.")
    with open(args.input, "r") as f:
        return json.load(f)


# --------------------------------------------------
def verify_dir(args: Args) -> None:
    """ Verify that the data directory and maps directory exist """    
    print(f"Validating data directory: {args.data}")
    if not args.data.exists():
        raise ValueError(f"Data directory '{args.data}' does not exist.")
    print(f"Validating maps directory: {args.data / 'maps'}")
    if not (args.data / "maps").exists():
        raise ValueError(f"Maps directory '{args.data / 'maps'}' does not exist.")


# --------------------------------------------------
def main() -> None:
    """ Run SAMap with command-line input """
    args = get_args()    
    
    verify_dir(args) # Ensure the data directory and maps directory exist
    data_dict = load_json(args) # Load the JSON file with species names and h5ad paths
    
    exit(0)  # Exit early if verification fails
    sm = SAMAP(
        load_json(args),
        f_maps = args.map_path,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
    )
    
    sm.run()
    save_samap(sm, 'example_data/samap_obj')


# --------------------------------------------------
if __name__ == '__main__':
    main()
