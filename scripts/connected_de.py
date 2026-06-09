#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-02-20
Version: 1.1.0
Purpose: Run Comprehensive Differential Expression analysis on a Pairwise species-species comparison
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
    input: Path #Path to the SAMap object
    pms: Path #Directory containing the pairwise mapping scores calculated upstream
    id1: str #id1 of pairwise comparison
    id2: str #id2 of pairwise comparison
    anno1: str #first annotation layer of pairwise comparison
    anno2: str #second annotation layer of pairwise comparison
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
        '-t', '--input',
        required=True,
        type=Path,
        help='Directory containing SAMap output object'
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
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    args = parser.parse_args()
    return Args(args.input, args.pms, args.id1, args.id2, args.anno1, args.anno2, args.output_dir)


# --------------------------------------------------
class ConnectedClusterDEAnalysis(object):
    def __init__(self, sm: object, keys: dict, pms_df, align_thr=0.2, grouping_thr=None, max_group_size=8, 
                 de_method='wilcoxon', min_cells=10):
        log(f"Starting init", "INFO")

        self.sm = sm
        self.keys = keys
        self.pms_df = pms_df
        self.align_thr = align_thr
        self.grouping_thr = grouping_thr if grouping_thr is not None else align_thr
        self.max_group_size = max_group_size
        self.de_method = de_method
        self.min_cells = min_cells
        
        log(f"Using pre-calculated PMS scores with threshold {align_thr}...", "INFO")
        
        # Step 1: Filter PMS scores by threshold
        self.high_alignment_pairs = self._filter_pms_by_threshold()
        
        # Step 2: Find connected clusters
        self.connected_clusters = self._find_connected_clusters()
        
        # Step 3: Prepare cell type data for DE analysis
        self.cluster_cell_data = self._prepare_cluster_cell_data()
        
        # Step 4: Generate biological names for each group
        self.group_names = self._generate_all_group_names()
        
        log(f"Found {len(self.connected_clusters)} connected alignment groups", "INFO")
        self._print_cluster_summary()
    
    def _filter_pms_by_threshold(self):
        log(f"Starting filter_pms_by_threshold", "INFO")
        
        high_align = self.pms_df[self.pms_df['pms_alignment_score'] >= self.align_thr].copy()
    
        cluster_cols = [col for col in self.pms_df.columns if 'cluster' in col]
        species1_col = cluster_cols[0]
        species2_col = cluster_cols[1]
        species1_prefix = species1_col.split('_')[0]
        species2_prefix = species2_col.split('_')[0]

        # Validate that the PMS file matches the species we're analyzing
        valid_prefixes = set(self.keys.keys())
        if not {species1_prefix, species2_prefix} == valid_prefixes:
            log(f"WARNING: PMS file species ({species1_prefix}, {species2_prefix}) do not match keys ({valid_prefixes})", "ERROR")
            raise ValueError(f"Wrong PMS file provided - contains {species1_prefix}/{species2_prefix} but expected {valid_prefixes}")

        high_align = high_align[
            (high_align[species1_col].apply(lambda x: species1_prefix in valid_prefixes)) &
            (high_align[species2_col].apply(lambda x: species2_prefix in valid_prefixes))
        ]

        log(f"Found {len(high_align)} cluster pairs above PMS threshold {self.align_thr}", "INFO")
        log(f"PMS score range for high-alignment pairs: {high_align['pms_alignment_score'].min():.4f} to {high_align['pms_alignment_score'].max():.4f}", "INFO")

        # Convert to the format expected by the rest of the code
        alignment_pairs = {}
        for _, row in high_align.iterrows():
            cluster1 = row[species1_col]
            cluster2 = row[species2_col]
            pms_score = row['pms_alignment_score']
            
            pair_key = f"{species1_prefix}_{cluster1};{species2_prefix}_{cluster2}"
            alignment_pairs[pair_key] = pms_score
        
        return alignment_pairs
    
    def _find_connected_clusters(self):
        log(f"Starting find_connected_clusters", "INFO")
        
        # Build graph from high-alignment pairs using grouping threshold
        graph = defaultdict(set)
        
        for pair_key, pms_score in self.high_alignment_pairs.items():
            if pms_score >= self.grouping_thr:
                ct1, ct2 = pair_key.split(';')
                graph[ct1].add(ct2)
                graph[ct2].add(ct1)
        
        log(f"Built graph with {len(graph)} nodes using grouping threshold {self.grouping_thr}", "INFO")
        
        # Find connected components
        visited = set()
        clusters = []
        
        def dfs(node, current_cluster):
            if node in visited or len(current_cluster) >= self.max_group_size:
                return
            visited.add(node)
            current_cluster.add(node)
            for neighbor in graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, current_cluster)
        
        for node in graph:
            if node not in visited:
                current_cluster = set()
                dfs(node, current_cluster)
                if len(current_cluster) > 1:  # Only keep multi-cell-type clusters
                    clusters.append(current_cluster)
        
        return clusters
    
    def _prepare_cluster_cell_data(self):
        log(f"Starting prepare_cluster_cell_data", "INFO")

        cluster_data = {}
        
        for i, cluster in enumerate(self.connected_clusters):
            cluster_id = f"cluster_{i}"
            cluster_data[cluster_id] = {
                'cell_types': cluster,
                'species_data': {},
                'pms_scores': {}
            }
            
            # Get PMS scores within this cluster
            cluster_pms = {}
            for ct1 in cluster:
                for ct2 in cluster:
                    if ct1 != ct2:
                        pair_key = f"{ct1};{ct2}"
                        reverse_key = f"{ct2};{ct1}"
                        if pair_key in self.high_alignment_pairs:
                            cluster_pms[pair_key] = self.high_alignment_pairs[pair_key]
                        elif reverse_key in self.high_alignment_pairs:
                            cluster_pms[reverse_key] = self.high_alignment_pairs[reverse_key]
            
            cluster_data[cluster_id]['pms_scores'] = cluster_pms
            
            # Group by species and fix cell type names
            species_groups = defaultdict(list)
            for cell_type in cluster:
                log(f"    Processing cell type: {cell_type}", "DEBUG")
                
                if '_' in cell_type:
                    parts = cell_type.split('_', 1)  # Split only on FIRST underscore
                    species_id = parts[0]
                    ct_part = parts[1] if len(parts) > 1 else parts[0]
                    
                    # Remove .0 suffix if present
                    if ct_part.endswith('.0'):
                        ct_num = ct_part[:-2]  # Remove '.0'
                    else:
                        ct_num = ct_part
                    
                    log(f"      Parsed as species={species_id}, cell_type={ct_num}", "DEBUG")
                    species_groups[species_id].append(ct_num)
            
            # Get cell data for each species in this cluster
            for species_id, cell_type_nums in species_groups.items():
                log(f"    Checking species {species_id} with cell types: {cell_type_nums}", "DEBUG")
                
                if species_id in self.keys:
                    adata = self.sm.sams[species_id].adata
                    ct_col = self.keys[species_id] 
                    
                    log(f"      Looking in column '{ct_col}'", "DEBUG")
                    log(f"      Available cell types in data: {sorted(adata.obs[ct_col].unique())}", "DEBUG")
                    
                    # Get all cells belonging to cell types in this cluster
                    cluster_mask = adata.obs[ct_col].isin(cell_type_nums)
                    n_cluster_cells = cluster_mask.sum()
                    
                    log(f"      Found {n_cluster_cells} cells matching cell types: {cell_type_nums}", "DEBUG")
                    
                    if n_cluster_cells >= self.min_cells:
                        cluster_data[cluster_id]['species_data'][species_id] = {
                            'adata': adata[cluster_mask].copy(),
                            'cell_types': cell_type_nums,
                            'n_cells': n_cluster_cells
                        }
                        log(f"      Added {species_id} to cluster with {n_cluster_cells} cells", "DEBUG")
                    else:
                        log(f"      Skipping {species_id} - insufficient cells ({n_cluster_cells} < {self.min_cells})", "DEBUG")
                else:
                    log(f"      Species {species_id} not found in SAMap object", "DEBUG")
        
        return cluster_data
    
    def _generate_all_group_names(self):
        """Generate biological names for all alignment groups"""
        log(f"Starting generate_all_group_names", "INFO")

        group_names = {}
        
        print("\nGenerating biological names for alignment groups...")
        for group_id in self.cluster_cell_data.keys():
            name = self._generate_consensus_name(group_id)
            group_names[group_id] = name
            log(f"  {group_id} → {name}", "INFO")
        
        return group_names
    
    def _generate_consensus_name(self, group_id):
        """Simple but always meaningful naming with fallback strategies"""
        log(f"Starting generate_consensus_name", "INFO")

        cluster_data = self.cluster_cell_data[group_id]
        cell_types = cluster_data['cell_types']
        
        parsed_names = []
        for ct in cell_types:
            if '_' in ct:
                type_name = '_'.join(ct.split('_')[1:]).replace('.0', '')
                if not type_name.isdigit():
                    parsed_names.append(type_name.lower())
        
        from collections import Counter
        if parsed_names:
            name_counts = Counter(parsed_names)
            most_common_name, count = name_counts.most_common(1)[0]
            
            if count > 1:
                # Clear consensus
                name_part = most_common_name
            elif len(parsed_names) <= 3:
                # Few types - list them all
                name_part = '+'.join(sorted(set(parsed_names)))
            else:
                # Many diverse types - use first + count
                name_part = f"{parsed_names[0]}_plus{len(set(parsed_names))-1}"
        else:
            name_part = "numeric_types"
        
        avg_pms = np.mean(list(cluster_data['pms_scores'].values())) if cluster_data['pms_scores'] else 0
        
        return f"{name_part}_pms{avg_pms:.2f}"
    
    def get_group_name(self, group_id):
        """Get the biological name for a group"""
        return self.group_names.get(group_id, group_id)
        
    def _print_cluster_summary(self):
        log(f"Starting print_cluster_summary", "INFO")

        log("\nAlignment Groups Summary:", "INFO")
        for cluster_id, data in self.cluster_cell_data.items():
            group_name = self.get_group_name(cluster_id)
            log(f"  {group_name} ({cluster_id}):", "INFO")
            log(f"    Cell types: {sorted(list(data['cell_types']))}", "INFO")
            
            # Show PMS scores within cluster
            if data['pms_scores']:
                pms_values = list(data['pms_scores'].values())
                avg_pms = np.mean(pms_values)
                log(f"    Average PMS score: {avg_pms:.4f}", "INFO")
                log(f"    PMS pairs in group:", "INFO")
                for pair, score in data['pms_scores'].items():
                    log(f"      {pair}: {score:.4f}", "INFO")
            
            for species_id, species_data in data['species_data'].items():
                log(f"    {species_id}: {len(species_data['cell_types'])} cell types, {species_data['n_cells']} cells", "INFO")
    
    def run_differential_expression_analysis(self, cluster_id, min_pct=0.5, min_logfc=1.0):

        log(f"Starting run_differential_expression_analysis", "INFO")

        if cluster_id not in self.cluster_cell_data:
            log(f"Group {cluster_id} not found!", "ERROR")
            return None
        
        cluster_data = self.cluster_cell_data[cluster_id]
        de_results = {}
        
        group_name = self.get_group_name(cluster_id)
        log(f"\nRunning Marker Gene analysis for {group_name} ({cluster_id})...", "INFO")
        
        for species_id, species_data in cluster_data['species_data'].items():
            log(f"  Analyzing {species_id}...", "INFO")
            
            adata = species_data['adata']
            
            # Get background cells (all other cell types, not just unclustered ones)
            full_adata = self.sm.sams[species_id].adata
            ct_col = self.keys[species_id]
            
            # Get only the current cluster's cell types for this species
            current_cluster_cts = set(species_data['cell_types'])
            
            # Background: ALL other cell types (including those in other clusters)
            all_cts = set(full_adata.obs[ct_col].astype(str).unique())
            background_cts = all_cts - current_cluster_cts
            
            log(f"    Current cluster cell types: {sorted(current_cluster_cts)}", "INFO")
            log(f"    Background cell types: {sorted(background_cts)}", "INFO")
            
            if len(background_cts) == 0:
                log(f"    No background cell types for {species_id}", "ERROR")
                continue
            
            background_mask = full_adata.obs[ct_col].astype(str).isin(background_cts)
            n_background_cells = background_mask.sum()
            
            if n_background_cells < self.min_cells:
                log(f"    Insufficient background cells for {species_id}: {n_background_cells}", "ERROR")
                continue
            
            # Prepare cluster and background cells
            cluster_cells = adata.copy()
            background_cells = full_adata[background_mask].copy()
            
            log(f"    Cluster cells: {cluster_cells.shape[0]}", "INFO")
            log(f"    Background cells: {background_cells.shape[0]}", "INFO")
            
            # Add group labels before concatenation
            cluster_cells.obs['in_cluster'] = 'True'
            background_cells.obs['in_cluster'] = 'False'
            
            # Also add a source column for batch correction if needed, In case we want to buff this function up in the future for timecourse DE, etc
            cluster_cells.obs['source'] = 'cluster'
            background_cells.obs['source'] = 'background'
            
            # Concatenate
            try:                
                combined_adata = cluster_cells.concatenate(
                    background_cells,
                    batch_key='batch',
                    batch_categories=['cluster', 'background'],
                    index_unique=None
                )
                
                log(f"    Combined data shape: {combined_adata.shape}", "INFO")
                log(f"    in_cluster values: {combined_adata.obs['in_cluster'].value_counts()}", "INFO")
                
                # Run Marker / DE analysis
                de_genes = self._run_de_scanpy(
                    combined_adata, 
                    groupby='in_cluster', 
                    reference='False',
                    min_pct=min_pct,
                    min_logfc=min_logfc
                )
                
                de_results[species_id] = {
                    'de_genes': de_genes,
                    'n_cells_cluster': len(cluster_cells),
                    'n_cells_background': len(background_cells),
                    'cell_types_in_cluster': species_data['cell_types'],
                    'background_cell_types': list(background_cts),
                    'pms_scores': cluster_data['pms_scores']
                }
                
                log(f"    Found {len(de_genes)} Marker genes in {species_id}")
                
            except Exception as e:
                log(f"    Marker analysis failed for {species_id}: {e}", "ERROR")
                import traceback
                traceback.print_exc()
                continue
        
        return de_results if de_results else None

    def _run_de_scanpy(self, adata, groupby, reference, min_pct=0.5, min_logfc=1.0):

        log(f"Starting run_de_scanpy", "INFO")
       
        # Ensure we have the right data
        if adata.raw is not None:
            adata_de = adata.raw.to_adata()
            adata_de.obs = adata.obs.copy()
        else:
            adata_de = adata.copy()
        
        # Ensure the groupby column is categorical
        if groupby in adata_de.obs.columns:
            adata_de.obs[groupby] = adata_de.obs[groupby].astype('category')
        else:
            raise ValueError(f"Column {groupby} not found in obs!")
        
        # Ensure we have sufficient cells in each group
        group_counts = adata_de.obs[groupby].value_counts()
        log(f"    Group counts: {dict(group_counts)}", "INFO")
        
        if len(group_counts) < 2:
            raise ValueError("Need at least 2 groups for DE analysis")
        
        if group_counts.min() < 5:
            raise ValueError(f"Insufficient cells in groups: {dict(group_counts)}")
        
        try:
            log(f"    Running scanpy rank_genes_groups...", "INFO")
            # Run differential expression - scanpy handles sparse matrices
            sc.tl.rank_genes_groups(
                adata_de,
                groupby=groupby,
                reference=reference,
                method=self.de_method,
                use_raw=False,
                key_added='de_analysis'
            )
            
            # Extract results
            if reference == 'False':
                group_name = 'True'
            else:
                group_name = reference
            
            result = sc.get.rank_genes_groups_df(adata_de, group=group_name, key='de_analysis')
            
            # DEBUG: Check what columns are actually available
            log(f"    Available columns: {list(result.columns)}", "INFO")
            
            # Filter by thresholds - need to handle different column names
            if len(result) > 0:
                # Handle different possible column names for log fold change
                logfc_col = None
                for col_name in ['logfoldchanges', 'lfc', 'log2FoldChange', 'logFC']:
                    if col_name in result.columns:
                        logfc_col = col_name
                        break
                
                # Handle different possible column names for percentage expressing
                pct_col = None
                for col_name in ['pct_nz_group', 'pct_expressing', 'pct', 'pct_nonzero']:
                    if col_name in result.columns:
                        pct_col = col_name
                        break
                
                # Apply filters based on available columns
                mask = pd.Series([True] * len(result))
                
                if logfc_col is not None:
                    mask = mask & (abs(result[logfc_col]) >= min_logfc)
                    log(f"    After logFC filter: {mask.sum()} genes", "INFO")
                
                if pct_col is not None:
                    mask = mask & (result[pct_col] >= min_pct)
                    log(f"    After pct filter: {mask.sum()} genes", "INFO")
                
                filtered_result = result[mask].copy()
                
                log(f"    Final filtered results: {len(filtered_result)} genes", "INFO")
                
                # Add FDR correction if we have p-values
                pval_col = None
                for col_name in ['pvals', 'pval', 'p_val', 'pvalue']:
                    if col_name in filtered_result.columns:
                        pval_col = col_name
                        break
                
                if pval_col is not None and len(filtered_result) > 0:
                    filtered_result['fdr'] = false_discovery_control(filtered_result[pval_col])

                # Filter by FDR threshold
                before_fdr = len(filtered_result)
                filtered_result = filtered_result[filtered_result['fdr'] <= 0.1].copy()
                log(f"    After FDR ≤ 0.1 filter: {len(filtered_result)} genes (removed {before_fdr - len(filtered_result)})", "INFO")
                
                # Sort by absolute log fold change if available
                if logfc_col is not None and len(filtered_result) > 0:
                    filtered_result = filtered_result.sort_values(logfc_col, key=abs, ascending=False)
                
                return filtered_result
            else:
                log("    No results returned from scanpy", "ERROR")
                return pd.DataFrame()
                
        except Exception as e:
            log(f"    Scanpy Marker analysis failed: {e}", "ERROR")
            import traceback
            traceback.print_exc()
            raise e
    
    def run_all_clusters_de_analysis(self, **de_params):

        log(f"Starting run_all_clusters_de_analysis", "INFO")

        import gc
        all_results = {}
        
        for cluster_id in self.cluster_cell_data.keys():
            results = self.run_differential_expression_analysis(cluster_id, **de_params)
            if results:
                all_results[cluster_id] = results
            
            gc.collect()
        
        return all_results
    
    def export_de_results(self, cluster_id, de_results, output_dir=None):
        """Export DE results with biological group name to group-specific folder"""

        log("starting export_de_results", "INFO") 

        group_name = self.get_group_name(cluster_id)
        
        # Create group-specific folder
        if output_dir is None:
            group_folder = group_name
        else:
            group_folder = os.path.join(output_dir, group_name)
        
        os.makedirs(group_folder, exist_ok=True)
        
        cluster_data = self.cluster_cell_data[cluster_id]
        
        # Export summary to group folder
        summary_path = os.path.join(group_folder, f"{group_name}_summary.txt")
        with open(summary_path, 'w') as f:
            f.write(f"Marker Gene Analysis: {group_name}\n")
            f.write(f"Internal ID: {cluster_id}\n")
            f.write(f"Cell types in group: {sorted(list(cluster_data['cell_types']))}\n")
            
            # Add PMS scores
            f.write(f"\nPMS Alignment Scores within group:\n")
            for pair, score in cluster_data['pms_scores'].items():
                f.write(f"  {pair}: {score:.4f}\n")
            
            f.write(f"\nSpecies Analysis Results:\n")
            for species_id, results in de_results.items():
                f.write(f"\n{species_id}:\n")
                f.write(f"  Cell types in group: {results['cell_types_in_cluster']}\n")
                f.write(f"  Background cell types: {results['background_cell_types'][:10]}...\n")
                f.write(f"  Group cells: {results['n_cells_cluster']}\n")
                f.write(f"  Background cells: {results['n_cells_background']}\n")
                f.write(f"  Marker genes: {len(results['de_genes'])}\n")
        
        # Export DE genes for each species to group folder
        for species_id, results in de_results.items():
            de_df = results['de_genes']
            de_path = os.path.join(group_folder, f"{group_name}_{species_id}_DE_genes.csv")
            de_df.to_csv(de_path, index=False)
        
        log(f"Exported DE results for {group_name} to {group_folder}/", "INFO")
        
    def _find_gene_name_column(self, de_genes_df):
        possible_cols = ['names', 'gene', 'gene_name', 'gene_id', 'symbol']
        
        for col in possible_cols:
            if col in de_genes_df.columns:
                return col
        
        # If no standard column found, use the first column
        log(f"Warning: No standard gene column found. Available columns: {list(de_genes_df.columns)}", "ERROR")
        return de_genes_df.columns[0]


    def export_barcode_alignment_labels(self, samap, output_path="barcode_alignment_families.csv"):
        """Export per-barcode alignment family labels for Loupe file creation."""
        log("Exporting barcode alignment family labels", "INFO")
        
        records = []
        for cluster_id, cluster_data in self.cluster_cell_data.items():
            group_name = self.get_group_name(cluster_id)
            for species_id, species_data in cluster_data['species_data'].items():
                adata = samap.sams[species_id].adata
                ct_col = self.keys[species_id]
                cell_types = species_data['cell_types']
                mask = adata.obs[ct_col].isin(cell_types)
                barcodes = adata.obs_names[mask]
                for barcode in barcodes:
                    records.append({'barcode': barcode, 'alignment_family': group_name})
        
        df = pd.DataFrame(records).set_index('barcode')
        df.to_csv(output_path)
        log(f"Exported {len(df)} barcode alignment labels to '{output_path}'", "INFO")
        return df

    def get_alignment_confusion_table(self):
        tables = {}
        
        for species_id, ct_col in self.keys.items():
            adata = self.sm.sams[species_id].adata
            all_cell_types = adata.obs[ct_col].unique()
            
            # Build mapping: cell_type → alignment_family
            ct_to_family = {ct: "Unassigned" for ct in all_cell_types}
            
            for cluster_id, cluster_data in self.cluster_cell_data.items():
                if species_id in cluster_data['species_data']:
                    group_name = self.get_group_name(cluster_id)
                    for ct in cluster_data['species_data'][species_id]['cell_types']:
                        ct_to_family[ct] = group_name
            
            # Count cells per (cell_type, family) pair
            obs = adata.obs[[ct_col]].copy()
            obs['alignment_family'] = obs[ct_col].map(ct_to_family)
            
            table = pd.crosstab(obs[ct_col], obs['alignment_family'])
            tables[species_id] = table
        
        return tables


# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the SAMap output object
    3. Loads the provided pairwise keys into the environment
    4. Runs and Exports a Differential Expression analysis for each Alignment Group within the pairiwise comparison
    """
    args = get_args()

    with open(args.input, "rb") as f:
        samap = pickle.load(f)
    
    # Get the keys you're using
    keys = {args.id1: args.anno1, args.id2: args.anno2}
    log(f"Pairiwise comparison is using keys dictionary: {keys}", "INFO")

    #Read in the pms table from the previous module
    pms_df = pd.read_csv(args.pms)

    # Initialize the analysis
    analysis = ConnectedClusterDEAnalysis(
        sm=samap,
        keys=keys,
        pms_df=pms_df,
        align_thr=0.3,
        grouping_thr=0.4,
        de_method='wilcoxon',
        min_cells=50
    )

    # Need this specific object in a downstream module
    with open('analysis.pkl', 'wb') as f:
        pickle.dump(analysis, f)

    analysis.export_barcode_alignment_labels(samap, "barcode_alignment_families.csv")

    # Calculate Celltype - Alignment Family Confusion Matrix
    tables = analysis.get_alignment_confusion_table()
    for species_id, table in tables.items():
        table.to_csv(f"{species_id}_alignment_confusion.csv")

    # Run DE analysis for all groups
    all_de_results = analysis.run_all_clusters_de_analysis(
        min_pct=0.5,
        min_logfc=3.0
    )

    # Need this specific object in a downstream module
    with open('all_de_results.pkl', 'wb') as f:
        pickle.dump(all_de_results, f)

    # Export DE results
    for cluster_id, de_results in all_de_results.items():
        analysis.export_de_results(
            cluster_id, 
            de_results, 
            output_dir="Grouping_Analysis" 
        )


    
# --------------------------------------------------
if __name__ == '__main__':
    main()


