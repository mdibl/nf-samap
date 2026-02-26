#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-02-18
Version: 1.1.0
Purpose: Provide Summary-Level Analysis on SAMap output object
"""

import argparse
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # For any plotting
import seaborn as sns  # For heatmaps in commented sections
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
import gc
import os
from collections import defaultdict
import scanpy as sc
from scipy import sparse
from scipy.stats import false_discovery_control
from log_utils import log



# --------------------------------------------------
def get_args() -> Args:
    """
    Parse and return command-line arguments.

    Returns:
        Args: A named tuple containing parsed command-line arguments for the output SAMap object, two id values, and two annotation layer values.
    """
    parser = argparse.ArgumentParser(
        description='Run Summary-Level Anlysis of the output SAMap object',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-t', '--input',
        required=True,
        type=Path,
        help='Directory containing SAMap output object'
    )

    parser.add_argument(
        '-i', '--id1',
        required=True,
        type=str,
        help='id1 of pairwise comparison'
    )

    parser.add_argument(
        '-d', '--id2',
        required=False,
        type=str,
        help='id2 of pairwise comparison'
    )

    parser.add_argument(
        '-a', '--anno1',
        required=True,
        type=Path,
        help='first annotation layer of pairwise comparison'
    )
    
    parser.add_argument(
        '-n', '--anno2',
        required=False,
        type=str,
        help='second annotation layer of pairwise comparison'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    args = parser.parse_args()
    return Args(args.input, args.id1, args.id2, args.anno1, args.anno2, args.output_dir)


# --------------------------------------------------
# Can I designate "keys" as a dict? Or do I have to loop through the function and pass id, annotation column individually?
def clean_annotations(keys: dict, samap: object) -> object:     
    for species_id, annotation_col in keys.items():
        
        # 1. Fix individual SAM
        sam = samap.sams[species_id]    
        sam.adata.obs[annotation_col] = (
            sam.adata.obs[annotation_col]
            .astype(str)
            .str.replace('_', '')
            .str.replace('/', '')
            .str.replace('-', '')
            .astype('category')
        )    
        # 2. Fix combined SAMap adata
        combined_col = f"{species_id}_{annotation_col.replace('.', '')}"
        
        if combined_col in samap.samap.adata.obs.columns:        
            samap.samap.adata.obs[combined_col] = (
                samap.samap.adata.obs[combined_col]
                .astype(str)
                .str.replace('_', '')
                .astype('category')
            )
        else:
            # Try without replacing '.' in the column name
            alt_combined_col = f"{species_id}_{annotation_col}"
            if alt_combined_col in samap.samap.adata.obs.columns:
                samap.samap.adata.obs[alt_combined_col] = (
                    samap.samap.adata.obs[alt_combined_col]
                    .astype(str)
                    .str.replace('_', '')
                    .astype('category')
                )
            else:
                log(f"\nWARNING: Neither {combined_col} nor {alt_combined_col} found!", "INFO")
                log(f"Available columns in combined adata:", "INFO")
                log(f"{[col for col in samap.samap.adata.obs.columns if species_id in col]}", "INFO")

    return(samap)


# --------------------------------------------------
def analyze_cluster_alignment(sm: object, keys, threshold=0.01, filter_threshold_pct=1.0, create_heatmaps=True) -> None:
    """
    Analyze cross-species cluster alignment from SAMap object.
    
    Parameters:
    -----------
    sm : SAMap object
        SAMap object containing mapping data
    keys : dict
        Dictionary mapping species prefixes to cluster column names
        e.g., {'mu': 'CellTypesIdent', 'ze': 'HarmonyRes.0.3'}
    threshold : float
        Threshold for cell-level alignment scores
    filter_threshold_pct : float
        Minimum percentage above threshold for "good" mappings
    create_heatmaps : bool
        Whether to create visualization heatmaps
    
    Returns:
    --------
    pd.DataFrame : Summary dataframe with alignment statistics
    """    
    # Get species and validate
    species_list = list(keys.keys())
    if len(species_list) != 2:
        raise ValueError(f"Expected 2 species, got {len(species_list)}")
    
    species1, species2 = species_list[0], species_list[1]
    
    # Create column names dynamically
    col1 = f"{species1}_{keys[species1].replace('.', '')}"
    col2 = f"{species2}_{keys[species2].replace('.', '')}"
    
    # Extract data
    combined_adata = sm.samap.adata #Need to specify which species?
    connectivity = combined_adata.obsp['connectivities']
    species_labels = combined_adata.obs['species'].values
    
    # Get species masks and indices
    mask1 = species_labels == species1
    mask2 = species_labels == species2
    indices1 = np.where(mask1)[0]
    indices2 = np.where(mask2)[0]
    
    # Extract cross-species mapping matrix
    cross_species_mapping = connectivity[indices1, :][:, indices2]
    
    # Get cluster information
    clusters1 = combined_adata.obs[col1][mask1].values
    clusters2 = combined_adata.obs[col2][mask2].values
    
    unique_clusters1 = sorted(np.unique(clusters1))
    unique_clusters2 = sorted(np.unique(clusters2))
    
    # Analyze cluster mappings
    log(f"\n=== Running Cluster Mapping Analysis ===", "INFO")
    log(f"Using threshold: {threshold}", "INFO")
    log(f"Analyzing {len(unique_clusters1)} {species1} clusters × {len(unique_clusters2)} {species2} clusters...", "INFO")
    
    summary_df = _compute_cluster_statistics(
        cross_species_mapping, clusters1, clusters2, 
        unique_clusters1, unique_clusters2, 
        species1, species2, threshold
    )
    # Create heatmaps
    if create_heatmaps:
        _create_alignment_heatmaps(summary_df, species1, species2)
    
    # Save results
    summary_df.to_csv('samap_alignment_analysis.csv', index=False)
        
    # Save filtered results
    good_mappings = summary_df[summary_df['pct_above_threshold'] > filter_threshold_pct]
    good_mappings.to_csv('samap_good_alignments.csv', index=False)
    
    log(f"Total results: {full_path}", "INFO")
    log(f"Filtered results (>{filter_threshold_pct}% above threshold): {filtered_path} ({len(good_mappings)} pairs)", "INFO")


def _compute_cluster_statistics(cross_species_mapping, clusters1, clusters2, 
                                unique_clusters1, unique_clusters2, 
                                species1, species2, threshold):
    """Compute alignment statistics for all cluster pairs."""
    
    results = []
    
    for c1 in unique_clusters1:
        cluster1_indices = np.where(clusters1 == c1)[0]
        
        for c2 in unique_clusters2:
            cluster2_indices = np.where(clusters2 == c2)[0]
            
            # Extract submatrix for this cluster pair
            submatrix = cross_species_mapping[np.ix_(cluster1_indices, cluster2_indices)]
            total_pairs = submatrix.shape[0] * submatrix.shape[1]
            
            # Calculate statistics
            if submatrix.nnz > 0:
                dense = submatrix.toarray()
                above_threshold = np.sum(dense > threshold)
                mean_score = np.mean(dense)
                median_score = np.median(dense)
                max_score = np.max(dense)
                std_score = np.std(dense)
            else:
                above_threshold = 0
                mean_score = median_score = max_score = std_score = 0.0
            
            results.append({
                f'{species1}_cluster': c1,
                f'{species2}_cluster': c2,
                f'{species1}_cells_in_cluster': len(cluster1_indices),
                f'{species2}_cells_in_cluster': len(cluster2_indices),
                'total_cell_pairs': total_pairs,
                'pairs_above_threshold': above_threshold,
                'pct_above_threshold': (above_threshold / total_pairs) * 100,
                'mean_alignment_score': mean_score,
                'median_alignment_score': median_score,
                'max_alignment_score': max_score,
                'std_alignment_score': std_score,
                'total_nonzero_pairs': submatrix.nnz
            })
    
    return pd.DataFrame(results)


def _create_alignment_heatmaps(summary_df, species1, species2) -> None:
    """Create visualization heatmaps for alignment metrics as separate files."""
    
    col1 = f'{species1}_cluster'
    col2 = f'{species2}_cluster'
    
    metrics = ['pct_above_threshold', 'mean_alignment_score', 'pairs_above_threshold']
    titles = ['Percentage Above Threshold', 'Mean Alignment Score', 'Number of Pairs Above Threshold']
    formats = ['.1f', '.3f', '.0f']
    filenames = ['pct_above_threshold', 'mean_alignment_score', 'pairs_above_threshold']
    
    for metric, title, fmt, filename in zip(metrics, titles, formats, filenames):
        # Create individual figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        pivot_df = summary_df.pivot(index=col1, columns=col2, values=metric)
        
        # Sort indices numerically if possible
        try:
            pivot_df = pivot_df.reindex(sorted(pivot_df.index, key=lambda x: float(x)))
            pivot_df = pivot_df.reindex(sorted(pivot_df.columns, key=lambda x: float(x)), axis=1)
        except:
            pass

        sns.heatmap(pivot_df, annot=True, fmt=fmt, ax=ax, 
                   cmap='viridis', cbar_kws={'shrink': 0.8})
        ax.set_title(title)
        ax.set_xlabel(f'{species2.upper()} Clusters')
        ax.set_ylabel(f'{species1.upper()} Clusters')
        
        plt.tight_layout()
        
        # Save if output directory provided
        plt.savefig(f'heatmap_{filename}.png', dpi=300, bbox_inches='tight')

        plt.close()

# --------------------------------------------------

# ENHANCE PMS
def extract_pms_scores(sm, keys, threshold=0.2):
    """
    Extract pairwise mapping scores (PMS) between species clusters from SAMap object.
    
    Parameters:
    -----------
    sm : SAMap object
        SAMap object containing mapping data
    keys : dict
        Dictionary mapping species prefixes to cluster column names
        e.g., {'mu': 'CellTypesIdent', 'ze': 'HarmonyRes.0.3'}
    threshold : float
        Threshold for cell-level alignment statistics
    project_id : str, optional
        Project identifier for output files
    
    Returns:
    --------
    tuple : (pms_df, enhanced_pms_df)
        Simple and enhanced PMS dataframes
    """
    # Get mapping scores
    hms, pms = get_mapping_scores(sm=sm, keys=keys)
        
    # Get combined AnnData
    combined_adata = sm.samap.adata
    species_labels = combined_adata.obs['species'].values
    
    # Convert keys to list for consistent ordering
    species_list = list(keys.keys())
    if len(species_list) != 2:
        raise ValueError(f"Expected 2 species, got {len(species_list)}")
    
    species1, species2 = species_list[0], species_list[1]
    
    # Create column names dynamically
    col1 = f"{species1}_{keys[species1].replace('.', '')}"
    col2 = f"{species2}_{keys[species2].replace('.', '')}"
    
    # Get species masks and clusters
    mask1 = species_labels == species1
    mask2 = species_labels == species2
    
    clusters1 = combined_adata.obs[col1][mask1].values
    clusters2 = combined_adata.obs[col2][mask2].values
    
    unique_clusters1 = sorted(np.unique(clusters1))
    unique_clusters2 = sorted(np.unique(clusters2))
    
    log(f"Species 1 ({species1}): {len(unique_clusters1)} clusters", "INFO")
    log(f"Species 2 ({species2}): {len(unique_clusters2)} clusters", "INFO")
    
    # FIXED: Extract cross-species PMS using actual index/column names
    pms_results = []
    
    for c1 in unique_clusters1:
        for c2 in unique_clusters2:
            # Construct the full cluster names as they appear in the matrix
            full_name1 = f"{species1}_{c1}"
            full_name2 = f"{species2}_{c2}"
            
            # Get the actual value from the matrix using the names
            if full_name1 in pms.index and full_name2 in pms.columns:
                score = pms.loc[full_name1, full_name2]
            elif full_name2 in pms.index and full_name1 in pms.columns:
                # Try reverse order
                score = pms.loc[full_name2, full_name1]
            else:
                log(f"Warning: Could not find {full_name1} vs {full_name2} in PMS matrix", "ERROR")
                score = 0.0
            
            pms_results.append({
                f'{species1}_cluster': c1,
                f'{species2}_cluster': c2,
                'pms_alignment_score': score
            })
    
    pms_df = pd.DataFrame(pms_results)
    
    log(f"\nTop 10 highest PMS alignment scores:", "INFO")
    top_pms = pms_df.nlargest(10, 'pms_alignment_score')
    for _, row in top_pms.iterrows():
        log(f"{species1}_{row[f'{species1}_cluster']} ↔ {species2}_{row[f'{species2}_cluster']}: PMS = {row['pms_alignment_score']:.6f}", "INFO")
    
    # Add cell-level statistics
    enhanced_pms_df = _add_cell_level_stats(
        pms_df, clusters1, clusters2, combined_adata, 
        species1, species2, mask1, mask2, threshold
    )
    
    # Save results
    pms_df.to_csv('pms_cluster_alignment_scores.csv', index=False)
    enhanced_pms_df.to_csv('pms_cluster_alignment_scores_enhanced.csv', index=False)
    
    # Summary statistics 
    log(f"\n=== PMS Summary Statistics ===", "INFO")
    log(f"Total cluster pairs: {len(pms_df)}")
    log(f"PMS score range: {pms_df['pms_alignment_score'].min():.6f} to {pms_df['pms_alignment_score'].max():.6f}", "INFO")
    log(f"Mean PMS score: {pms_df['pms_alignment_score'].mean():.6f}", "INFO")
    log(f"Median PMS score: {pms_df['pms_alignment_score'].median():.6f}", "INFO")
    log(f"Standard deviation: {pms_df['pms_alignment_score'].std():.6f}", "INFO")
    
    return pms_df, enhanced_pms_df


def _add_cell_level_stats(pms_df, clusters1, clusters2, combined_adata, 
                          species1, species2, mask1, mask2, threshold):
    """
    Add cell-level alignment statistics to PMS dataframe.
    
    Internal helper function for extract_pms_scores.
    """
    connectivity = combined_adata.obsp['connectivities']
    indices1 = np.where(mask1)[0]
    indices2 = np.where(mask2)[0]
    cross_species_mapping = connectivity[indices1, :][:, indices2]
    
    col1_name = f'{species1}_cluster'
    col2_name = f'{species2}_cluster'
    
    enhanced_results = []
    
    for _, row in pms_df.iterrows():
        cluster1 = row[col1_name]
        cluster2 = row[col2_name]
        
        # Get cell indices for this cluster pair
        cluster1_indices = np.where(clusters1 == cluster1)[0]
        cluster2_indices = np.where(clusters2 == cluster2)[0]
        
        # Calculate statistics
        if len(cluster1_indices) > 0 and len(cluster2_indices) > 0:
            submatrix = cross_species_mapping[np.ix_(cluster1_indices, cluster2_indices)]
            submatrix_dense = submatrix.toarray() if hasattr(submatrix, 'toarray') and submatrix.nnz > 0 else (
                submatrix if not hasattr(submatrix, 'toarray') else np.zeros(submatrix.shape)
            )
            
            if submatrix_dense.size > 0:
                mean_score = np.mean(submatrix_dense)
                median_score = np.median(submatrix_dense)
                std_score = np.std(submatrix_dense)
                pct_above = (np.sum(submatrix_dense > threshold) / submatrix_dense.size) * 100
            else:
                mean_score = median_score = std_score = pct_above = 0.0
        else:
            mean_score = median_score = std_score = pct_above = 0.0
        
        enhanced_results.append({
            col1_name: cluster1,
            col2_name: cluster2,
            'pms_alignment_score': row['pms_alignment_score'],
            'mean_cell_score': mean_score,
            'median_cell_score': median_score,
            'std_cell_score': std_score,
            'pct_above_threshold': pct_above
        })
    
    return pd.DataFrame(enhanced_results)





# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the SAMap output object
    3. Clean the SAMap output object
    4. Creates heatmaps and runs summary-level analysis
    """
    args = get_args()

    with open(args.input, "rb") as f:
        samap = pickle.load(f)
    
    # Get the keys you're using
    keys = {args.id1: args.anno1, args.id2: args.anno2}
    log(f"Pairiwise comparison is using keys dictionary: {keys}", "INFO")


    log("Fixing any potential issues in the Annotation Labels", "INFO")
    samap = clean_annotations(keys, samap)

    analyze_cluster_alignment(samap, keys, threshold=0.01, 
                              filter_threshold_pct=1.0, create_heatmaps=True)


    log("Finding Gene Pairs for Pairwise Species Comparison", "INFO")
    gpf = GenePairFinder(samap, keys=keys)
    gene_pairs = gpf.find_all(align_thr=0.3)
    gene_pairs.to_csv(f'GenePairs.csv', index=False)

    
    log("Performing PMS dataframe enchancement", "INFO")
    pms_df, enhanced_pms_df = extract_pms_scores(samap, keys, threshold=0.2)

    pms_df.to_csv('pms_tidy.csv', index=False)

    log("Writing Updated SAMap Object out to pickle format", "INFO")
    with open('samap_results_Cleaned.pkl', 'wb') as f:
        pickle.dump(samap, f)
# --------------------------------------------------
if __name__ == '__main__':
    main()


