#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP with command-line input
"""

import argparse
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple


class Args(NamedTuple):
    """ Command-line arguments """
    species1: str
    h5ad1: Path
    species2: str
    h5ad2: Path
    map_path: str
    


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Run SAMAP with command-line input',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('species1',
                        metavar='s1',
                        type=str,
                        help='First species name')
    
    parser.add_argument('h5ad1',
                        metavar='h5ad1',
                        type=Path,
                        help='First species h5ad file')
    
    parser.add_argument('species2',
                        metavar='s2',
                        type=str,
                        help='Second species name')    

    parser.add_argument('h5ad2',
                        metavar='h5ad2',
                        type=Path,
                        help='Second species h5ad file')
    
    parser.add_argument('--map-path',
                        metavar='PATH',
                        type=str,
                        default='example_data/maps/',
                        help='Path to orthology map directory')
    

    args = parser.parse_args()

    return Args(args.species1,
                args.h5ad1,
                args.species2,
                args.h5ad2,
                args.map_path)


# --------------------------------------------------
def main() -> None:
    """ Run SAMap with command-line input """
    args = get_args()
    
    sm = SAMAP(
        build_h5ad_dict(args),
        f_maps = args.map_path,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
    )
    
    sm.run()
    save_samap(sm, 'example_data/samap_output')

# --------------------------------------------------
def build_h5ad_dict(args: Args) -> dict:
    """ Build a dictionary of h5ad files for SAMAP """
    return {
        args.species1: str(args.h5ad1),
        args.species2: str(args.h5ad2),
    }
    

# --------------------------------------------------
if __name__ == '__main__':
    main()
