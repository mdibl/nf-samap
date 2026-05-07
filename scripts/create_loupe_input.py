#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-05-07
Version: 2.0.3
Purpose: Extract SAMap embedding and combine with raw counts for Loupe file creation
"""

import argparse
import pickle
import pandas as pd
import numpy as np
import anndata as ad
from typing import NamedTuple
from pathlib import Path
from log_utils import log


class Args(NamedTuple):
    """Command-line arguments for the script"""
    samap:         Path    # Path to the SAMap pickle file
    h5ads:         list    # Paths to the raw AnnData h5ad files, one per species
    species:       list    # List of species IDs to subset from SAMap object
    celltype_cols: list    # List of obs column names for cell type labels per species
    output:        Path    # Path to the output combined h5ad file

# --------------------------------------------------

def get_args() -> Args:
    parser = argparse.ArgumentParser(
        description='Extract SAMap embedding and combine with raw counts for Loupe file creation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--samap',
        required=True,
        type=Path,
        help='Path to the SAMap pickle file'
    )
    parser.add_argument(
        '--h5ads',
        required=True,
        type=Path,
        nargs='+',
        help='Paths to raw h5ad files, one per species in same order as --species')

    parser.add_argument(
        '--species',
        required=True,
        type=str,
        nargs='+',
        help='Species IDs (e.g. ax ze)'
    )
    parser.add_argument(
        '--celltype-cols',
        required=True,
        type=str,
        nargs='+',
        help='obs column names for cell type labels, one per species in same order as --species'
    )

    args = parser.parse_args()
    return Args(args.samap, args.h5ads, args.species, args.celltype_cols)

# --------------------------------------------------

def load_samap(samap_path: Path):
    """Load and return the SAMap object from a pickle file."""
    log(f" Loading SAMap object from '{samap_path}'", "INFO")
    with open(samap_path, 'rb') as f:
        sm = pickle.load(f)
    log(f" Loaded SAMap object with species: {sm.ids}", "INFO")
    return sm

# --------------------------------------------------


def subset_samap(sm, species: list):
    """Subset the combined SAMap AnnData to the specified species."""
    log(f" Subsetting SAMap AnnData to species: {species}", "INFO")
    adata = sm.samap.adata
    adata_subset = adata[adata.obs['species'].isin(species)].copy()
    log(f" Subsetted to {adata_subset.n_obs} cells", "INFO")
    for sp in species:
        count = (adata_subset.obs['species'] == sp).sum()
        log(f"   {sp}: {count} cells", "INFO")
    return adata_subset

# --------------------------------------------------

def build_cell_type_labels(adata_subset, species: list, celltype_cols: list) -> pd.Series:
    """Build prefixed cell type labels for each species."""
    log("Building prefixed cell type labels", "INFO")
    labels = pd.Series(index=adata_subset.obs_names, dtype=str)
    for sp, col in zip(species, celltype_cols):
        mask = adata_subset.obs['species'] == sp
        labels[mask] = sp + '_' + adata_subset.obs.loc[mask, col].astype(str)
        log(f"   {sp}: using column '{col}'", "INFO")
    return labels

# --------------------------------------------------

def load_and_combine_h5ads(h5ad_paths: list, species: list, target_barcodes: pd.Index) -> ad.AnnData:
    """Load raw h5ad files, combine with outer join, and subset to target barcodes."""
    adatas = []
    for path, sp in zip(h5ad_paths, species):
        log(f" Loading h5ad for '{sp}' from '{path}'", "INFO")
        a = ad.read_h5ad(path)
        log(f"   {sp}: {a.n_obs} cells x {a.n_vars} genes", "INFO")
        adatas.append(a)

    log("Concatenating h5ad files with outer join", "INFO")
    combined = ad.concat(adatas, axis=0, join='outer', fill_value=0)
    log(f" Combined: {combined.n_obs} cells x {combined.n_vars} genes", "INFO")

    missing = set(target_barcodes) - set(combined.obs_names)
    if missing:
        raise ValueError(f"{len(missing)} SAMap barcodes not found in combined h5ad. First 5: {list(missing)[:5]}")

    log(f"Subsetting combined h5ad to SAMap cell order", "INFO")
    return combined[target_barcodes].copy()

# --------------------------------------------------

def main() -> None:
    args = get_args()

    sm = load_samap(args.samap)
    adata_subset = subset_samap(sm, args.species)

    adata_subset.obs['cell_type_labeled'] = build_cell_type_labels(
        adata_subset, args.species, args.celltype_cols
    )

    combined = load_and_combine_h5ads(args.h5ads, args.species, adata_subset.obs_names)

    log(f"Attaching SAMap UMAP and metadata", "INFO")
    combined.obsm['X_umap']              = adata_subset.obsm['X_umap']
    combined.obs['cell_type_labeled']    = adata_subset.obs['cell_type_labeled'].values
    combined.obs['species']              = adata_subset.obs['species'].values

    combined.write_h5ad(f"{args.species[0]}_{args.species[1]}_LoupeInput.h5ad")

# --------------------------------------------------

if __name__ == '__main__':
    main()