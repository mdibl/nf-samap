#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-09
Purpose: Visualize SAMAP results from pickle
"""

import argparse
import os
import pickle
from typing import NamedTuple, Optional
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import matplotlib.pyplot as plt

class Args(NamedTuple):
    input: str
    output_dir: Optional[str]


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Visualize SAMAP results from pickle',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar='PKL',
                        required=True,
                        help='Input SAMAP pickle file')

    parser.add_argument('-o', '--output-dir',
                        metavar='DIR',
                        default='output/',
                        help='Optional output directory for visualizations')

    args = parser.parse_args()

    return Args(args.input, args.output_dir)


# --------------------------------------------------
def main() -> None:
    """ Visualize SAMAP results from pickle"""

    args = get_args()

    # Ensure output directory exists
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"Output directory ensured: {args.output_dir}")

    # Load the pickle file
    with open(args.input, 'rb') as f:
        samap_obj = pickle.load(f)
    print(f"Loaded SAMAP object of type: {type(samap_obj)}")

    keys = {'pl':'cluster','hy':'Cluster','sc':'tissue'}
    highest_mapping_scores, pairwise_mapping_scores = get_mapping_scores(samap_obj, keys, n_top=0)

    # Save mapping dataframes to CSV
    hms_outfile = os.path.join(args.output_dir, "highest_mapping_scores.csv")
    highest_mapping_scores.to_csv(hms_outfile)
    print(f"Saved highest mapping scores to: {hms_outfile}")
    pms_outfile = os.path.join(args.output_dir, "pairwise_mapping_scores.csv")
    pairwise_mapping_scores.to_csv(pms_outfile)
    print(f"Saved pairwise mapping scores to: {pms_outfile}")
    
    # Scatter plot
    samap_obj.scatter()
    sc_outfile = os.path.join(args.output_dir, "samap_scatter.png")
    plt.savefig(sc_outfile, dpi=300)
    print(f"Saved scatter plot to {sc_outfile}")
    plt.close()

# --------------------------------------------------
if __name__ == '__main__':
    main()
