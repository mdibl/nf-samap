#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2025-06-16
Version: 1.0.0
Purpose: Use the .csv and sparse Matrix files generated in the previous preprocessing module to create an AnnData object
"""

import os
import sys

# CRITICAL: Set performance environment variables BEFORE any imports
os.environ['OMP_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'
os.environ['MKL_NUM_THREADS'] = '4'
os.environ['NUMEXPR_MAX_THREADS'] = '4'
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['MPLBACKEND'] = 'Agg'

# Only disable numba caching, NOT JIT compilation (for performance)
os.environ['NUMBA_DISABLE_CACHING'] = '1'

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
    counts: Path        # Path to that species' Sparse Matrix containing counts information
    obs: Path           # Path to the .obs file created for that species
    feats: Path         # Path to the .feats file created for that species
    output_dir: Path    # Path to the output directory


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
        '--counts',
        required=True,
        type=Path,
        help='Path to the Sparse Matrix containing counts information for that species'
    )

    parser.add_argument(
        '--obs',
        required=True,
        type=Path,
        help='Path to the .obs file created for that species'
    )

    parser.add_argument(
        '--feats',
        required=True,
        type=Path,
        help='Path to the .feats file created for that species'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    args = parser.parse_args()
    return Args(args.id, args.counts, args.obs, args.feats, args.output_dir)


# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the necessary .csv and Sparse Matrix files.
    3. Creates the initial AnnData Object.
    4. Runs standard preprocessing on the AnnData object to prepare it for SAMap.
    5. Saves the h5ad object.
    """

    log("Loading arguments", "INFO")
    args = get_args()

    # 1. Load matrix (cells x genes, from .mtx file)
    log("Attempting to load the Sparse Counts Matrix", "INFO")
    X = mmread(args.counts).tocsr()
    log("Successfully loaded the Sparse Counts Matrix!", "INFO")
    log(f"Dimensions of Counts: {X.shape[0]} vs {X.shape[1]}", "INFO")

    # 2. Load metadata

    #Loading Cell IDs + metadata
    log("Attempting to load the Cell IDs, desired Annotation Layer, and other relevant cell metadata", "INFO")
    obs = pd.read_csv(args.obs, header=None)
    log("Successfully loaded!", "INFO")


    obs = obs.rename(columns=obs.iloc[0]).drop(0).reset_index(drop=True)
    obs = obs.set_index("cell_id")  # set cell IDs as index, IMPORTANT DO NOT REMOVE, SAMap doesn't know what to call the cells otherwise


    #Loading Gene IDs
    log("Attempting to load the Gene IDs", "INFO")
    var = pd.read_csv(args.feats, header=None)
    log("Successfully loaded!", "INFO")


    var.columns = ["gene"]
    

    # 3. Check alignment

    #Checking Cells vs Cell IDs
    log("Checking to see if there exists a discrepancy in number of Cell IDs vs. number of Cells", "INFO")
    if X.shape[0] != obs.shape[0]:
        log(f"Cells mismatch: {X.shape[0]} Cells vs {obs.shape[0]} Cell IDs", "ERROR")

    #Checking Genes vs Gene IDs
    log("Checking to see if there exists a discrepancy in number of Gene IDs vs. number of Genes", "INFO")
    if X.shape[1] == var.shape[0]:
        log(f"Genes mismatch: {X.shape[1]} Genes vs {var.shape[0]} Gene IDs", "ERROR")

    # 4. Build AnnData
    log("Attempting to Initialize the AnnData object", "INFO")
    var = var.set_index("gene")
    adata = ad.AnnData(X=X, obs=obs, var=var)
    log("Successfully Initialized Anndata Object!", "INFO")


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
        n_genes=3000,
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