#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-02-28
Version: 1.1.0
Purpose: Use upstream analyses to draw >2 species conclusions about alignment
"""

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # For any plotting
import gc
import os
from collections import defaultdict
import scanpy as sc
from scipy import sparse
from scipy.stats import false_discovery_control
from log_utils import log
from typing import NamedTuple
from pathlib import Path
import argparse



# --------------------------------------------------

class Args(NamedTuple):
    analysis: Path #Path to the analysis object generated upstream
    output_dir: Path #Path to the output directory


def get_args() -> Args:
    """
    Parse and return command-line arguments.

    Returns:
        Args: A named tuple containing parsed command-line arguments for the output SAMap object, two id values, and two annotation layer values.
    """
    parser = argparse.ArgumentParser(
        description='Run Summary-Level Anlysis of the output SAMap object',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_arguments(
        '-a', '--analysis',
        required=True,
        type=Path,
        nargs="+",
        help='Path to the output directory'
    )

    
    parser.add_argument(
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    args = parser.parse_args()
    return Args(args.output_dir)


# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads all Analysis Folders from previous downstream analyses as a single list
    3. Uses all pair-linked analyses to draw >2 species conclusions
    """
    args = get_args()

    keys = {args.id1: args.anno1, args.id2: args.anno2} #CHANGE TO >2 SPECIES, needs to be flexible

    #HERE'S WHAT I WANT TO DO IN THIS MODULE:
    # 1.) n>2-species pms map + alignment families
    # 2.) Gene Triangles, Gene 3-way Horizontal comparisons (nuance here in exact relationship, 
    # also present 1a -> 2 -> 3 -> 1b as putative homolog if "1a <-/-> 1b" as a gene pair)
    #   a. Also integrate Marker gene information into gene triangle / gene horizontal comparison table.
    #   b. Make filter for # of diff species existing in >2 species gene relationship
    #
    

# --------------------------------------------------
if __name__ == '__main__':
    main()


