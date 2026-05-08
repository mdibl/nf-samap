#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2025-06-16
Version: 1.0.0
Purpose: Create SAM objects from the AnnData input data, preprocess them using SAM's preprocess functions
"""

import os
import sys

# CRITICAL: Set performance environment variables BEFORE any imports
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['MPLBACKEND'] = 'Agg'

# Only disable numba caching, NOT JIT compilation (for performance)
# os.environ['NUMBA_DISABLE_CACHING'] = '1'

from log_utils import log
log("Loaded Log_utils", "INFO")

log("Loading pandas...", "INFO")
import pandas as pd

log("Loading scipy...", "INFO")
import scipy.sparse as sp
from scipy.io import mmread

log("Loading anndata...", "INFO")
import anndata as ad
log("Loaded anndata", "INFO")

log("Loading scanpy...", "INFO")  # This is the critical one
import scanpy as sc
log("Loaded scanpy", "INFO")

log("Loading samalg...", "INFO")
import samalg  # make sure samalg is installed (this is the SAM library)
log("Loaded samalg", "INFO")

log("Loading remaining packages...", "INFO")
import argparse
from pathlib import Path
log("Loaded pathlib", "INFO")
from typing import NamedTuple

log("ALL IMPORTS SUCCESSFUL!", "INFO")



log("Loaded all Packages!", "INFO")


class Args(NamedTuple):
    """ Command-line arguments for the script"""
    
    id: str             # Species ID for the sample being processed
    anndata: Path       # Path to the AnnData Object containing the expression data, either computed directly upstream or provided as input sample
    vargenes: int         # Number of variables genes the user would like to be considered in AnnData Preprocessing (default = 3000)

# --------------------------------------------------
def get_args() -> Args:
    """
    Parse and return command-line arguments.

    Returns:
        Args: A named tuple containing parsed command-line arguments for id, counts, obs, feats, and output_dir.
    """
    parser = argparse.ArgumentParser(
        description=' Use the .csv and sparse Matrix files generated in the previous preprocessing module to create an AnnData object',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--id',
        required=True,
        type=str,
        help='Species ID for the sample being processed'
    )

    parser.add_argument(
        '--anndata',
        required=True,
        type=Path,
        help='Path to the AnnData object for that species'
    )
    parser.add_argument(
        '--vargenes',
        required=True,
        type=int,
        help='# of variable genes to be considered'
    )


    args = parser.parse_args()
    return Args(args.id, args.anndata, args.vargenes)

    # --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the necessary .csv and Sparse Matrix files.
    3. Creates the initial SAM Object.
    4. Runs standard preprocessing on the AnnData object to prepare it for SAMap.
    """

    log("Loading arguments", "INFO")
    args = get_args()

    adata = args.anndata
    var_genes = args.vargenes

    # 5. Wrap AnnData in SAM object
    sam = samalg.SAM(adata)

    # 6. Preprocess (adjust params to match test h5ad)
    log("Attempting to run preprocessing on the AnnData Object", "INFO")
    sam.preprocess_data()

    # 7. Run SAM (these parameters match the test run_args you shared)
    sam.run(
        k=20,
        distance="cosine",
        projection="umap",
        npcs=150,
        n_genes=var_genes,
        max_iter=10,
        seed=None,
        sparse_pca=False,
        weight_PCs=False,
        weight_mode="combined",
        verbose=True
    )
    log("Successfully preprocessed the AnnData Object!", "INFO")

    # 8. Save as h5ad (AnnData v0.7.8 compatible)
    log("Saving the Preprocessed AnnData object", "INFO")
    sam.adata.write(f"{args.output_dir}/{args.id}_preprocessed.h5ad")
    log("Successfully saved!", "INFO")

# --------------------------------------------------
if __name__ == "__main__":
    main()