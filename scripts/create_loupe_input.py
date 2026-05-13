#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-05-07
Version: 2.0.5
Purpose: Extract SAMap embedding and combine with raw counts for Loupe file creation
"""

import argparse
import pickle
import pandas as pd
import numpy as np
import anndata as ad
from scipy.io import mmwrite
from typing import NamedTuple
from pathlib import Path
from log_utils import log


class Args(NamedTuple):
    """Command-line arguments for the script"""
    samap:         Path         # Path to the SAMap pickle file
    h5ads:         list         # Paths to the raw AnnData h5ad files, one per species
    species:       list         # List of species IDs to subset from SAMap object
    celltype_cols: list         # List of obs column names for cell type labels per species
    alignment_families: Path    # Path to barcode alignment family labels CSV


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
        help='Paths to raw h5ad files, one per species in same order as --species'
    )
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
    parser.add_argument(
        '--alignment-families',
        required=False,
        type=Path,
        default=None,
        help='Path to barcode alignment family labels CSV'
    )

    args = parser.parse_args()

    if len(args.h5ads) != len(args.species):
        raise ValueError(f"Number of --h5ads ({len(args.h5ads)}) must match number of --species ({len(args.species)})")
    if len(args.celltype_cols) != len(args.species):
        raise ValueError(f"Number of --celltype-cols ({len(args.celltype_cols)}) must match number of --species ({len(args.species)})")

    return Args(args.samap, args.h5ads, args.species, args.celltype_cols, args.alignment_families)

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
        prefixed_col = f"{sp}_{col}"
        actual_col = prefixed_col if prefixed_col in adata_subset.obs.columns else col
        if actual_col not in adata_subset.obs.columns:
            raise KeyError(
                f"Column '{col}' not found for species '{sp}'. "
                f"Tried '{prefixed_col}' and '{col}'. "
                f"Available columns: {adata_subset.obs.columns.tolist()}"
            )
        mask = adata_subset.obs['species'] == sp
        labels[mask] = sp + '_' + adata_subset.obs.loc[mask, actual_col].astype(str)
        log(f"   {sp}: using column '{actual_col}'", "INFO")
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

    log("Subsetting combined h5ad to SAMap cell order", "INFO")
    return combined[target_barcodes].copy()

# --------------------------------------------------

def main() -> None:
    args = get_args()
    prefix = f"{args.species[0]}_{args.species[1]}"

    sm           = load_samap(args.samap)
    adata_subset = subset_samap(sm, args.species)

    adata_subset.obs['cell_type_labeled'] = build_cell_type_labels(
        adata_subset, args.species, args.celltype_cols
    )

    combined = load_and_combine_h5ads(args.h5ads, args.species, adata_subset.obs_names)

    log("Attaching SAMap UMAP and metadata", "INFO")
    combined.obsm['X_umap']           = adata_subset.obsm['X_umap']
    combined.obs['cell_type_labeled'] = adata_subset.obs['cell_type_labeled'].values
    combined.obs['species']           = adata_subset.obs['species'].values

    # Merge alignment families into metadata if provided
    if args.alignment_families:
        log("Loading alignment family labels", "INFO")
        af_df = pd.read_csv(args.alignment_families, index_col='barcode')
        combined.obs['alignment_family'] = af_df.reindex(combined.obs_names)['alignment_family'].fillna('Unassigned')
        log(f"Assigned alignment families to {(combined.obs['alignment_family'] != 'Unassigned').sum()} cells", "INFO")
    else:
        combined.obs['alignment_family'] = 'Unassigned'

    # Export count matrix in Matrix Market format (transposed: genes x cells for R)
    log("Exporting count matrix in Matrix Market format", "INFO")
    mmwrite(f"{prefix}_counts.mtx", combined.X.T)

    # Export barcodes and features
    log("Exporting barcodes and features", "INFO")
    pd.Series(combined.obs_names).to_csv(f"{prefix}_barcodes.csv", index=False, header=False)
    pd.Series(combined.var_names).to_csv(f"{prefix}_features.csv", index=False, header=False)

    # Export UMAP
    log("Exporting UMAP coordinates", "INFO")
    umap_df = pd.DataFrame(
        combined.obsm['X_umap'],
        index=combined.obs_names,
        columns=['UMAP_1', 'UMAP_2']
    )
    umap_df.to_csv(f"{prefix}_umap.csv")

    # Export metadata
    log("Exporting metadata", "INFO")
    combined.obs[['cell_type_labeled', 'species', 'alignment_family']].to_csv(f"{prefix}_meta.csv")

    log(f" Done. All files written with prefix '{prefix}'", "INFO")

# --------------------------------------------------

if __name__ == '__main__':
    main()