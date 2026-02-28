#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-02-20
Version: 1.1.0
Purpose: Use upstream pms and Differential Expression calculations to compare gene pair expression levels + more
"""

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # For any plotting
import seaborn as sns  # For heatmaps in commented sections
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles,
                            get_mapping_scores)
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
    genepairs: Path #Path to pairwise gene pairs calculated upstream
    diff: Path #Path to directory of DE data calculated upstream
    pms: Path #Directory containing the pairwise mapping scores calculated upstream
    id1: str #id1 of pairwise comparison
    id2: str #id2 of pairwise comparison
    anno1: str #first annotation layer of pairwise comparison
    anno2: str #second annotation layer of pairwise comparison
    analysis: Path #Path to the analysis object generated upstream
    all_de_results: Path #Path to the all_de_results object generated upstream
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

    parser.add_argument(
        '-g', '--genepairs',
        required=True,
        type=Path,
        help='Path to pairwise gene pairs calculated upstream'
    )

    parser.add_argument(
        '-f', '--diff',
        required=True,
        type=Path,
        help='Path to directory of DE data calculated upstream'
    )


    parser.add_argument(
        '-p', '--pms',
        required=True,
        type=Path,
        help='Directory containing the pairwise mapping scores calculated upstream'
    )

    parser.add_argument(
        '-i', '--id1',
        required=True,
        type=str,
        help='id1 of pairwise comparison'
    )

    parser.add_argument(
        '-d', '--id2',
        required=True,
        type=str,
        help='id2 of pairwise comparison'
    )

    parser.add_argument(
        '-a', '--anno1',
        required=True,
        type=str,
        help='first annotation layer of pairwise comparison'
    )
    
    parser.add_argument(
        '-n', '--anno2',
        required=True,
        type=str,
        help='second annotation layer of pairwise comparison'
    )
    
    parser.add_argument(
        '-y', '--analysis',
        required=True,
        type=Path,
        help='Path to the analysis object generated upstream'
    )
    parser.add_argument(
        '-s', '--all_de_results',
        required=True,
        type=Path,
        help='Path to the all_de_results object generated upstream'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    args = parser.parse_args()
    return Args(args.genepairs, args.diff, args.pms, args.id1, args.id2, args.anno1, args.anno2, args.analysis, args.all_de_results, args.output_dir)


# --------------------------------------------------

def add_cluster_info_to_dataframes(pairs_sep, connected_clusters, analysis):
    """
    Add cluster info to dataframes with biological group names
    
    Parameters:
    -----------
    pairs_sep : list
        List of dataframes with gene pairings
    connected_clusters : list
        List of connected clusters from analysis
    analysis : ConnectedClusterDEAnalysis
        Analysis object with group names
    
    Returns:
    --------
    list : Dataframes with 'group' column containing biological names
    """
    filtered_pairs_with_clusters = []
    
    for df in pairs_sep:
        # Get the first two column names (these are the cell type pairs with species prefix)
        col1 = df.columns[0]  # e.g., "mu_Excitatory_neurons" 
        col2 = df.columns[1]  # e.g., "ze_Excitatory_neurons"
        
        # Extract just the cell type names (everything after first underscore)
        if '_' in col1:
            parts1 = col1.split('_', 1)  # Split ONLY on first underscore
            species1 = parts1[0]
            firstCol = parts1[1] if len(parts1) > 1 else col1
        else:
            firstCol = col1
            species1 = col1
            
        if '_' in col2:
            parts2 = col2.split('_', 1)  # Split ONLY on first underscore
            species2 = parts2[0]
            secondCol = parts2[1] if len(parts2) > 1 else col2
        else:
            secondCol = col2
            species2 = col2
        
        # Remove .0 if present
        firstCol = firstCol.replace('.0', '')
        secondCol = secondCol.replace('.0', '')
        
        # Check which cluster this pair belongs to
        found_cluster = None
        cluster_cell_types = None
        group_name = None
        
        for cluster_idx, cluster in enumerate(connected_clusters):
            cluster_nums = set()
            for ct in cluster:
                if '_' in ct:
                    # Extract cell type name (everything after first underscore)
                    ct_parts = ct.split('_', 1)
                    ct_num = ct_parts[1] if len(ct_parts) > 1 else ct_parts[0]
                    ct_num = ct_num.replace('.0', '')
                    cluster_nums.add(ct_num)
            
            # Check if both numbers are in this cluster
            if firstCol in cluster_nums and secondCol in cluster_nums:
                found_cluster = f"cluster_{cluster_idx}"
                cluster_cell_types = sorted(list(cluster_nums))
                group_name = analysis.get_group_name(found_cluster)
                break
        
        # If we found a cluster, add the dataframe with cluster info
        if found_cluster is not None:
            # Make a copy of the dataframe
            df_with_cluster = df.copy()
            
            # Add group column with biological name
            df_with_cluster['group'] = group_name
            
            filtered_pairs_with_clusters.append(df_with_cluster)
            
            log(f"Added pair ({firstCol}, {secondCol}) to {group_name}", "INFO")
        else:
            log(f"WARNING: No cluster found for pair ({firstCol}, {secondCol})", "ERROR")
            for cluster_idx, cluster in enumerate(connected_clusters):
                cluster_nums = set()
                for ct in cluster:
                    if '_' in ct:
                        ct_parts = ct.split('_', 1)
                        ct_num = ct_parts[1] if len(ct_parts) > 1 else ct_parts[0]
                        ct_num = ct_num.replace('.0', '')
                        cluster_nums.add(ct_num)

    
    return filtered_pairs_with_clusters
    

# --------------------------------------------------

def add_de_info_to_dataframes(filtered_pairs_sep, all_de_results, keys, outdir, analysis):
    """
    Add DE information to dataframes and save with biological group names
    
    Parameters:
    -----------
    filtered_pairs_sep : list
        List of dataframes with group information
    all_de_results : dict
        DE results keyed by cluster_id (e.g., 'cluster_0')
    keys : dict
        Species keys dictionary
    outdir : str
        Output directory path
    analysis : ConnectedClusterDEAnalysis
        Analysis object with group names
    
    Returns:
    --------
    list : Enhanced dataframes with DE information
    """
    enhanced_dataframes = []
    
    species_keys = list(keys.keys()) if isinstance(keys, dict) else keys
    log(f"Number of input dataframes: {len(filtered_pairs_sep)}", "INFO")

    for df in filtered_pairs_sep:
        group_name = df['group'].iloc[0]
        
        # Find the corresponding cluster_id from the group name
        cluster_id = None
        for cid, gname in analysis.group_names.items():
            if gname == group_name:
                cluster_id = cid
                break
        
        if cluster_id is None or cluster_id not in all_de_results:
            log(f"Warning: Group {group_name} not found in Marker results", "ERROR")
            enhanced_dataframes.append(df.copy())
            continue
        
        enhanced_df = df.copy()
        
        # Initialize new columns
        for species_key in species_keys:
            enhanced_df[f'{species_key}_logFC'] = pd.NA
            enhanced_df[f'{species_key}_pval_adj'] = pd.NA
        
        cluster_de_results = all_de_results[cluster_id]
        log(f"\n=== Processing {group_name} ({cluster_id}) ===", "INFO")
        
        # Get the gene columns from the pairing dataframe
        gene_columns = enhanced_df.columns[:2]
        
        # Process each row
        matches_found = {species_key: 0 for species_key in species_keys}
        
        for idx, row in enhanced_df.iterrows():
            # Match each gene column to its species
            for gene_col in gene_columns:
                gene_name = str(row[gene_col])
                
                # Determine which species this gene belongs to by checking the prefix
                matched_species = None
                for species_key in species_keys:
                    if gene_name.startswith(f"{species_key}_"):
                        matched_species = species_key
                        break
                
                # If we identified the species and have DE results for it
                if matched_species and matched_species in cluster_de_results:
                    species_de_genes = cluster_de_results[matched_species]['de_genes']
                    
                    if len(species_de_genes) > 0:
                        gene_match = species_de_genes[species_de_genes['names'] == gene_name]
                        
                        if not gene_match.empty:
                            enhanced_df.at[idx, f'{matched_species}_logFC'] = gene_match['logfoldchanges'].iloc[0]
                            enhanced_df.at[idx, f'{matched_species}_pval_adj'] = gene_match['pvals_adj'].iloc[0]
                            matches_found[matched_species] += 1
        
        # Construct the output directory path using biological name
        save_dir = f"{outdir}/{group_name}/{group_name}_GenePairMarkers"
        os.makedirs(save_dir, exist_ok=True)
        
        # Create filename from first two column names
        filename = enhanced_df.columns[0] + "-" + enhanced_df.columns[1]
        
        # Save the file
        enhanced_df.to_csv(f"{save_dir}/{filename}.GenePairMarkers.csv", index=False)
        
        enhanced_dataframes.append(enhanced_df)
        log(f"Saved to {save_dir}/{filename}.GenePairMarkers.csv", "INFO")
    
    return enhanced_dataframes


# --------------------------------------------------

def combine_dfs(enhanced_pairs, keys):
    """
    Combine enhanced dataframes into a single stacked dataframe
    
    Parameters:
    -----------
    enhanced_pairs : list
        List of enhanced dataframes
    keys : dict
        Species keys dictionary
    
    Returns:
    --------
    pd.DataFrame : Combined dataframe with standardized columns
    """
    # Reset all column names to generic placeholders
    enhanced_pairs = [df.set_axis(range(df.shape[1]), axis=1) for df in enhanced_pairs]

    stacked_df = pd.concat(enhanced_pairs, ignore_index=True)
    
    # Dynamically create column names based on keys
    species_keys = list(keys.keys())
    column_names = [
        species_keys[0],  # First species
        species_keys[1],  # Second species
        "pairing_pval1",
        "pairing_pval2",
        "group",
        f"{species_keys[0]}_logFC",
        f"{species_keys[0]}_pval_adj",
        f"{species_keys[1]}_logFC",
        f"{species_keys[1]}_pval_adj"
    ]
    
    stacked_df.columns = column_names

    first_col = stacked_df.columns[0]
    stacked_df = stacked_df[~stacked_df[first_col].astype(str).str.lower().isin(['nan', 'none', ''])]
    return stacked_df


# --------------------------------------------------

def compress_dfs(combined_df, output_dir, keys, analysis): #Need analysis in this module?
    """
    Compress dataframes and save with biological group names
    
    Parameters:
    -----------
    combined_df : pd.DataFrame
        Combined dataframe with all gene pairs
    output_dir : str
        Output directory path
    keys : dict
        Species keys dictionary
    analysis : ConnectedClusterDEAnalysis
        Analysis object with group names
    
    Returns:
    --------
    dict : Dictionary of compressed dataframes keyed by biological group name
    """
    dfs_by_group = {value: subset for value, subset in combined_df.groupby('group')}
    compressed_dict = {}
    
    species_keys = list(keys.keys())
    keys_list = [keys[sp] for sp in species_keys]  # Get the full key names
    
    for group_name, df in dfs_by_group.items():
        df = df.copy()
        
        # Coerce all columns except the first two to numeric
        cols_to_convert = df.columns[2:]
        df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors='coerce')
        
        df['duplicates'] = 1
        df['comb'] = df[species_keys[0]] + ";" + df[species_keys[1]].astype(str)

        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if 'duplicates' in numeric_cols:
            numeric_cols.remove('duplicates')
        if 'comb' in numeric_cols:
            numeric_cols.remove('comb')

        agg_dict = {col: 'mean' for col in numeric_cols}
        agg_dict['duplicates'] = 'sum'

        compressed = df.groupby('comb', as_index=False).agg(agg_dict)
        
        # Split the 'comb' column back into species columns
        split_cols = compressed['comb'].str.split(';', expand=True)
        compressed[species_keys[0]] = split_cols[0]
        compressed[species_keys[1]] = split_cols[1]
        
        # Reorder columns - PUT SPECIES COLUMNS FIRST
        other_cols = [col for col in compressed.columns if col not in species_keys + ['comb']]
        compressed = compressed[species_keys + other_cols]

        compressed_dict[group_name] = compressed

        filename = f"{group_name}/{group_name}_GenePairMarkers/{group_name}_compressed.csv"
        filepath = os.path.join(output_dir, filename)
        compressed.to_csv(filepath, index=False)

    return compressed_dict


# --------------------------------------------------

def find_unpaired_de_genes(all_de_results, compressed_dfs, analysis):
    """
    Find DE genes that don't have cross-species pairings
    Uses biological group names for tracking
    
    Parameters:
    -----------
    all_de_results : dict
        DE results keyed by cluster_id
    compressed_dfs : dict
        Compressed dataframes keyed by biological group name
    analysis : ConnectedClusterDEAnalysis
        Analysis object
    
    Returns:
    --------
    dict : Unpaired results keyed by biological group name
    """
    
    # Get species IDs dynamically from keys
    species_ids = list(analysis.keys.keys())    
    unpaired_results = {}
    
    for cluster_id, cluster_de_results in all_de_results.items():
        group_name = analysis.get_group_name(cluster_id)
        
        # Get paired genes dynamically for all species
        if group_name in compressed_dfs:
            pairing_df = compressed_dfs[group_name]
            paired_genes = {}
                            
            # Dynamically get paired genes for each species column
            for species_id in species_ids:
                if species_id in pairing_df.columns:
                    paired_genes[species_id] = set(pairing_df[species_id].dropna().astype(str).str.upper())
                else:
                    paired_genes[species_id] = set()
                    log(f"  {species_id} column NOT FOUND in pairing df", "ERROR")
        else:
            paired_genes = {species_id: set() for species_id in species_ids}
            log(f"{group_name}: No compressed df found", "ERROR")
        
        cluster_results = {}
        
        # Process each species in DE results
        for species_id, species_results in cluster_de_results.items():
            de_genes_df = species_results['de_genes']
            
            if len(de_genes_df) > 0:
                gene_col = analysis._find_gene_name_column(de_genes_df)
                de_gene_names = set(de_genes_df[gene_col].astype(str).str.upper())
                
                # Find unpaired genes
                unpaired_genes = de_gene_names - paired_genes.get(species_id, set())
                
                # Filter Marker dataframe to unpaired genes only
                unpaired_mask = de_genes_df[gene_col].astype(str).str.upper().isin(unpaired_genes)
                unpaired_df = de_genes_df[unpaired_mask].copy()
                
                cluster_results[species_id] = {
                    'unpaired_de_genes': unpaired_df,
                    'unpaired_count': len(unpaired_genes),
                    'total_de_count': len(de_gene_names)
                }
                
                log(f"{group_name} {species_id}: {len(unpaired_genes)}/{len(de_gene_names)} unpaired Marker genes", "INFO")
        
        unpaired_results[group_name] = cluster_results
    
    return unpaired_results

# -------------------------------------------------- Maybe Need bc of custom class in the pkl file being read in?

class ConnectedClusterDEAnalysis(object):
    def get_group_name(self, group_id):
        """Get the biological name for a group"""
        return self.group_names.get(group_id, group_id)
    
    def _find_gene_name_column(self, de_genes_df):
        possible_cols = ['names', 'gene', 'gene_name', 'gene_id', 'symbol']
        
        for col in possible_cols:
            if col in de_genes_df.columns:
                return col
        
        # If no standard column found, use the first column
        log(f"Warning: No standard gene column found. Available columns: {list(de_genes_df.columns)}", "ERROR")
        return de_genes_df.columns[0]
        
# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the path-linked arguments as objects
    3. Loads the provided pairwise keys into the environment
    4. Uses both the DE data and pre-discovered pairwise gene pairs to generate insights into gene pair expression
    """
    args = get_args()

    keys = {args.id1: args.anno1, args.id2: args.anno2}
    log(f"Pairwise comparison is using keys dictionary: {keys}", "INFO")

    with open(args.analysis, 'rb') as f: #Okay let's see if this works instead of deconvoluting what's in Grouping_Analysis
        analysis = pickle.load(f)
    log(f"Successfully loaded analysis object", "INFO")

    with open(args.all_de_results, 'rb') as f:  #Okay let's see if this works instead of deconvoluting what's in Grouping_Analysis
        all_de_results = pickle.load(f)
    log(f"Successfully loaded all_de_results object", "INFO")

    # -----------------------------------------------------------

    # Load and prepare gene pairs; reformat to make it easier to work with
    pairs = pd.read_csv(args.genepairs)
    pairs_list = [pairs.iloc[:, i:i+3] for i in range(0, len(pairs.columns), 3)]
    pairs_sep = []
    for group_df in pairs_list:
        if group_df.iloc[:, 0].isna().all():
                continue

        # Split first column
        colnames = group_df.columns[0].split(';')
        split_cols = group_df.iloc[:, 0].astype(str).str.split(';', expand=True)
        if split_cols.shape[1] != len(colnames):
            log(f"Warning: Column split mismatch for {group_df.columns[0]}", "ERROR")
            log(f"  Expected {len(colnames)} columns, got {split_cols.shape[1]}", "ERROR")
            continue
        split_cols.columns = colnames
        
        # Combine with remaining columns
        new_df = pd.concat([split_cols, group_df.iloc[:, 1:]], axis=1)
        pairs_sep.append(new_df)

    # -----------------------------------------------------------
    log(f"Attempting to add cluster info to Gene Pairs", "INFO")

    # Add cluster info with biological names
    filtered_pairs_sep = add_cluster_info_to_dataframes(
        pairs_sep, 
        analysis.connected_clusters, 
        analysis 
    )
    log(f"SUCCESSFULLY added cluster info to Gene Pairs", "INFO")

    log(f"Attempting to add DE information to Gene Pairs", "INFO")

    # Add DE information
    enhanced_pairs = add_de_info_to_dataframes(
        filtered_pairs_sep, 
        all_de_results, 
        keys, 
        "GroupingAnalysis",
        analysis 
    )

    log(f"SUCCESSFULLY added DE information to Gene Pairs", "INFO")


    #--------------- COMBINED AND COMPRESSED DATA ---------------
    log(f"Attempting to combine Gene Pairs into a single Dataframe", "INFO")

    combined_dfs = combine_dfs(enhanced_pairs, keys)  

    log(f"SUCCESSFULLY combined Gene Pairs into a single Dataframe", "INFO")

    log(f"Attempting to compress Gene Pairs df to have unique rows, count duplicates", "INFO")

    compressed_dfs = compress_dfs(
        combined_dfs, 
        "GroupingAnalysis", 
        keys,      
        analysis
    )
    log(f"SUCCESSFULLY compressed Gene Pairs df to have unique rows, count duplicates", "INFO")

    log("Keys in compressed_dfs:", list(compressed_dfs.keys()), "INFO")
    log("Sample key:", list(compressed_dfs.keys())[0] if compressed_dfs else "None", "INFO")

    #--------------- UNPAIRED DE GENE DISCOVERY ---------------

    log(f"Attempting to create a df of non-Gene Pair Marker Genes", "INFO")

    unpaired_de_results = find_unpaired_de_genes(all_de_results, compressed_dfs, analysis)

    log(f"SUCCESSFULLY created a df of non-Gene Pair Marker Genes", "INFO")

    # Export unpaired DE genes to group-specific folders
    for group_name, group_data in unpaired_de_results.items():
        group_folder = f"GroupingAnalysis/{group_name}/Unpaired_Marker"
        os.makedirs(group_folder, exist_ok=True)
        
        for species_id, species_data in group_data.items():
            if len(species_data['unpaired_de_genes']) > 0:
                filename = f"{group_folder}/unpaired_Marker_{species_id}.csv"
                species_data['unpaired_de_genes'].to_csv(filename, index=False)

    log("\nAnalysis complete!", "INFO")

# --------------------------------------------------
if __name__ == '__main__':
    main()


