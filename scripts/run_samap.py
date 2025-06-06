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
from typing import NamedTuple, TextIO


class Args(NamedTuple):
    """ Command-line arguments """
    species1: str
    h5ad1: Path
    map1to2: Path
    species2: str
    h5ad2: Path
    map2to1: Path


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
    
    parser.add_argument('map1to2',
                        metavar='map1to2',
                        type=Path,
                        help='Mapping from species 1 to species 2')
    
    parser.add_argument('species2',
                        metavar='s2',
                        type=str,
                        help='Second species name')    

    parser.add_argument('h5ad2',
                        metavar='h5ad2',
                        type=Path,
                        help='Second species h5ad file')
    
    parser.add_argument('map2to1', 
                        metavar='map2to1',
                        type=Path,
                        help='Mapping from species 2 to species 1')

    args = parser.parse_args()

    return Args(args.species1,
                args.h5ad1,
                args.map1to2,
                args.species2,
                args.h5ad2,
                args.map2to1)


# --------------------------------------------------
def main() -> None:
    """ Run SAMap with command-line input """

    args = get_args()
    
    
    species1 = args.species1
    h5ad1 = str(args.h5ad1)
    map1to2 = str(args.map1to2)
    
    species2 = args.species2
    h5ad2 = str(args.h5ad2)
    map2to1 = str(args.map2to1)
    

    print(f'Running SAMAP for {species1} vs {species2}')
    print(f'Input files: {h5ad1},              {h5ad2}')
    print(f'Mapping files: {map1to2},                {map2to1}')
    
    
    
    fn1 = h5ad1
    fn2 = h5ad2
    fn3 = 'example_data/hydra.h5ad'
    
    filenames = {'pl':fn1,'sc':fn2,'hy':fn3}
    sm = SAMAP(
        filenames,
        f_maps = 'example_data/maps/',
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
    )
    
    sm.run(pairwise=True)
    
    
    # keys = {'pl':'cluster','hy':'Cluster','sc':'tissue'}
    # D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
    save_samap(sm, 'example_data/samap_output')
    
    
    
    
# --------------------------------------------------
def run_samap(species1: str, h5ad1: TextIO, map1to2: TextIO,
            species2: str, h5ad2: TextIO, map2to1: TextIO) -> None:
    """ Run SAMAP with the provided arguments """
    
    # Pack the input parameters into a dictionary
    filenames = build_input_dict(species1, h5ad1, species2, h5ad2)
    
    # Create the SAMAP object
    sm = SAMAP(
        filenames,
        f_maps='/home/ryan/git/github/ryansonder/my-samap/example_data/maps/hypl'
    )
    
    sm.run(pairwise=True)
    samap = sm.samap
    
    
    # Run the SAMAP mapping
    
    
# --------------------------------------------------
    
def build_input_dict(species1: str, h5ad1: TextIO,
               species2: str, h5ad2: TextIO) -> dict:
    """ Build a dictionary of input parameters for SAMAP """
    return {
        species1: h5ad1,
        species2: h5ad2,
    }


# --------------------------------------------------
if __name__ == '__main__':
    main()
