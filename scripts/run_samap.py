#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP with command-line input
"""

import argparse
import json
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple


class Args(NamedTuple):
    """ Command-line arguments """
    map_path: str
    json_path: Path


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Run SAMAP with command-line input',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--map-path',
                        metavar='PATH',
                        type=str,
                        default='example_data/maps/',
                        help='Path to orthology map directory')

    parser.add_argument('--json-path',
                        type=Path,
                        required=True,
                        help='Path to JSON with species names and h5ad paths')

    args = parser.parse_args()

    return Args(args.map_path,
                args.json_path,)


# --------------------------------------------------
def main() -> None:
    """ Run SAMap with command-line input """
    args = get_args()    
    
    sm = SAMAP(
        load_json(args),
        f_maps = args.map_path,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
    )
    
    sm.run()
    save_samap(sm, 'example_data/samap_output')


# --------------------------------------------------
# def build_h5ad_dict(args: Args) -> dict:
#     """ Build a dictionary of h5ad files for SAMAP """
#     return {
#         args.species1: str(args.h5ad1),
#         args.species2: str(args.h5ad2),
#     }


# --------------------------------------------------
def load_json(args: Args) -> dict[str: Path]:
    with open(args.json_path, "r") as f:
        return json.load(f)

# --------------------------------------------------
if __name__ == '__main__':
    main()
