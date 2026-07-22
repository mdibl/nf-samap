#!/usr/bin/env python3
"""
Author : Markus Sujansky
Date   : 2026-07-22
Version: 0.1.0
Purpose: Per-cell-type divergence of BLAST-linked (homology) gene components.

For each cross-species cell-type pairing with PMS >= --pms-thr, walk each
connected component of the BLAST homology graph (sm.gnnm, thresholded on raw
bit score) and label every member gene -- on its OWN species' side of the
pairing -- as one of {up, down, flat, off}. A component is RETURNED for that
pairing iff its members do not all share a single label (internal incoherence).

Outputs (long / tidy; variable component size lives in rows, not columns):
  <id1>_<id2>_homology_divergence.csv       one row per (component x pairing x gene)
  <id1>_<id2>_pairing_provenance.csv        one row per (pairing x species side)
"""

import argparse
import pickle
from pathlib import Path
from typing import NamedTuple
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components

try:
    from log_utils import log
except Exception:  # allow standalone runs outside the pipeline container
    def log(msg, level="INFO"):
        print(f"[{level}] {msg}", flush=True)


# --------------------------------------------------
class Args(NamedTuple):
    """Command-line arguments for the script"""
    input:      Path    # Path to the SAMap pickle file (cleaned SAMap object)
    pms:        Path    # Path to pms_cluster_alignment_scores.csv from SUMMARY_SAMAP
    id1:        str     # First species ID of the pairwise comparison
    id2:        str     # Second species ID of the pairwise comparison
    anno1:      str     # obs column holding cell type labels for id1
    anno2:      str     # obs column holding cell type labels for id2
    output_dir: Path    # Directory to write the divergence and provenance CSVs
    pms_thr:    float   # Minimum pms_alignment_score for a cell-type pairing to be considered
    blast_thr:  float   # Minimum BLAST bit score (sm.gnnm edge weight) to keep a homology edge
    off_frac:   float   # Fraction-of-cells-expressing gate; below this a gene is labelled "off"
    fdr_thr:    float   # FDR (pvals_adj) threshold for up/down significance
    min_cells:  int     # Minimum cells in a foreground cell type required to run DE
    de_method:  str     # scanpy rank_genes_groups method (e.g. "wilcoxon")


def get_args() -> Args:

    parser = argparse.ArgumentParser(
        description="Per-cell-type divergence of BLAST-linked homology components",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-t',
        '--input',
        required=True,
        type=Path,
        help='SAMap output object (pickle)'
    )

    parser.add_argument(
        '-p',
        '--pms',
        required=True,
        type=Path,
        help='pms_cluster_alignment_scores.csv from SUMMARY_SAMAP'
    )

    parser.add_argument(
        '-i', 
        '--id1', 
        required=True, 
        type=str, 
        help='id1 of pairwise comparison'
    )

    parser.add_argument(
        '-d', 
        '--id2', 
        required=True, 
        type=str, 
        help='id2 of pairwise comparison'
    )

    parser.add_argument(
        '-a', 
        '--anno1', 
        required=True, 
        type=str, 
        help='annotation (obs col) for id1'
    )

    parser.add_argument(
        '-n', 
        '--anno2', 
        required=True, 
        type=str, 
        help='annotation (obs col) for id2'
    )

    parser.add_argument(
        '-o', 
        '--output_dir', 
        type=Path,
        default=Path('.')
    )

    parser.add_argument(
        '--pms-thr', 
        type=float, 
        default=0.3,
        help='minimum pms_alignment_score for a pairing to be considered'
    )

    parser.add_argument(
        '--blast-thr', 
        type=float, 
        default=0.0,
        help='minimum BLAST bit score (sm.gnnm edge weight) to keep an edge'
    )

    parser.add_argument(
        '--off-frac', 
        type=float, 
        default=0.10,
        help='fraction-of-cells-expressing gate; below this a gene is "off"'
    )

    parser.add_argument(
        '--fdr-thr', 
        type=float, 
        default=0.1,
        help='FDR (pvals_adj) threshold for up/down significance'
    )

    parser.add_argument(
        '--min-cells', 
        type=int, 
        default=50,
        help='minimum cells in a foreground cell type to run DE'
    )

    parser.add_argument('--de-method', type=str, default='wilcoxon')
    a = parser.parse_args()

    return Args(a.input, a.pms, a.id1, a.id2, a.anno1, a.anno2, a.output_dir,
                a.pms_thr, a.blast_thr, a.off_frac, a.fdr_thr, a.min_cells, a.de_method)


# --------------------------------------------------
def norm_ct(x) -> str:
    """Normalise a cell-type identifier: str, and drop a trailing '.0' float suffix.
    Mirrors the handling in connected_de.py so PMS cluster ids match adata.obs values."""
    s = str(x)
    return s[:-2] if s.endswith('.0') else s


def build_homology_components(sm, id1, id2, blast_thr):
    """Connected components of the BLAST gene graph, restricted to id1/id2 genes.

    sm.gnnm is a symmetric sparse gene x gene matrix of (reciprocal-filtered)
    BLAST bit scores, aligned to sm.gns. Only cross-species edges exist in it,
    so any surviving component spans both species by construction.

    Returns: list of components, each a list of full gene ids (e.g. 'ax_MMP13').
    """
    gns = np.asarray(sm.gns)
    log(f"sm.gns holds {gns.size} genes; first 5 look like: {list(gns[:5])}", "INFO")
    species_of = np.array([g.split('_', 1)[0] for g in gns])
    prefixes, counts = np.unique(species_of, return_counts=True)
    log(f"Gene prefixes present in sm.gns: {dict(zip(prefixes.tolist(), counts.tolist()))}", "INFO")
    keep = np.isin(species_of, [id1, id2])
    n1 = int((species_of == id1).sum())
    n2 = int((species_of == id2).sum())
    log(f"Genes matched by prefix -> {id1}: {n1}, {id2}: {n2}", "INFO")
    if n1 == 0 or n2 == 0:
        log(f"WARNING: expected prefixes '{id1}_' and '{id2}_' but one matched 0 genes. "
            f"Gene ids may not be prefixed as assumed (see sample above).", "ERROR")
    if keep.sum() == 0:
        raise ValueError(f"No genes with prefix {id1}/{id2} found in sm.gns")

    sub = sm.gnnm[keep][:, keep].tocsr()
    sub_gns = gns[keep]

    # threshold on bit score -> boolean adjacency
    adj = sub.copy()
    adj.data = (adj.data > blast_thr).astype(np.int8)
    adj.eliminate_zeros()

    n_comp, labels = connected_components(csgraph=adj, directed=False)
    groups = defaultdict(list)
    for idx, lab in enumerate(labels):
        groups[lab].append(sub_gns[idx])

    # keep only multi-gene components (a single gene cannot diverge)
    comps = [g for g in groups.values() if len(g) >= 2]
    if comps:
        sizes = np.array([len(c) for c in comps])
        log(f"Homology graph: {keep.sum()} {id1}/{id2} genes, {len(comps)} components "
            f">=2 members (bit-score > {blast_thr}); size min/median/max = "
            f"{sizes.min()}/{int(np.median(sizes))}/{sizes.max()}", "INFO")
    else:
        log(f"WARNING: 0 multi-gene components after thresholding at bit-score > {blast_thr}. "
            f"If this is unexpected, the threshold may be too high or sm.gnnm may be empty.", "ERROR")
    return comps


def qualifying_pairings(pms_df, id1, id2, pms_thr):
    """Return list of dicts: {celltype_1, celltype_2, pms} for pms >= threshold.
    celltype_1 is the id1 side, celltype_2 the id2 side (normalised)."""
    log(f"PMS file columns: {list(pms_df.columns)}", "INFO")
    cluster_cols = [c for c in pms_df.columns if 'cluster' in c]
    if len(cluster_cols) != 2:
        raise ValueError(f"Expected two *cluster* columns in PMS file, got {cluster_cols}")
    # figure out which cluster col is id1 vs id2 by prefix in the column name
    try:
        col1 = next(c for c in cluster_cols if c.split('_')[0] == id1)
        col2 = next(c for c in cluster_cols if c.split('_')[0] == id2)
    except StopIteration:
        raise ValueError(
            f"Could not match cluster columns {cluster_cols} to ids ({id1}, {id2}) by prefix. "
            f"Column naming may differ from the assumed '{{id}}_cluster'.")
    log(f"Matched PMS columns -> {id1}: '{col1}', {id2}: '{col2}'", "INFO")
    if 'pms_alignment_score' not in pms_df.columns:
        raise ValueError("PMS file has no 'pms_alignment_score' column")
    log(f"pms_alignment_score range: {pms_df['pms_alignment_score'].min():.4f} "
        f"to {pms_df['pms_alignment_score'].max():.4f} over {len(pms_df)} rows", "INFO")

    hi = pms_df[pms_df['pms_alignment_score'] >= pms_thr]
    pairings = []
    for _, r in hi.iterrows():
        pairings.append({
            'celltype_1': norm_ct(r[col1]),
            'celltype_2': norm_ct(r[col2]),
            'pms': float(r['pms_alignment_score']),
        })
    log(f"{len(pairings)} pairings with PMS >= {pms_thr}", "INFO")
    return pairings


def per_celltype_de(sm, species_id, ct_col, cell_type, de_method, min_cells):
    """One-vs-rest DE for a SINGLE cell type within one species.
    Returns (gene_records, provenance) or (None, provenance) if skipped.

    gene_records: dict gene_id -> {'lfc','fdr','frac','mean'} keyed WITHOUT species
                  prefix stripped (uses adata.var_names as-is; caller prefixes).
    """
    full = sm.sams[species_id].adata
    if ct_col not in full.obs.columns:
        raise ValueError(f"annotation column '{ct_col}' not in {species_id}.obs")

    obs_ct = full.obs[ct_col].astype(str).map(norm_ct)
    fg_mask = (obs_ct == cell_type).values
    n_fg = int(fg_mask.sum())
    n_bg = int((~fg_mask).sum())
    prov = {'n_cells': n_fg, 'n_background_cells': n_bg, 'median_depth': np.nan}

    # n_fg == 0 almost always means a name mismatch between the PMS cluster id
    # and adata.obs values -- surface both sides so it is diagnosable.
    if n_fg == 0:
        avail = sorted(obs_ct.unique().tolist())
        log(f"  [{species_id}] cell type '{cell_type}' matched 0 cells in obs['{ct_col}']. "
            f"Available (normalised) values: {avail[:20]}"
            f"{' ...' if len(avail) > 20 else ''}", "ERROR")
        return None, prov

    if n_fg < min_cells or n_bg < min_cells:
        log(f"  [{species_id}] skip cell type '{cell_type}': "
            f"fg={n_fg}, bg={n_bg} (min_cells={min_cells})", "INFO")
        return None, prov

    # use the same expression source connected_de.py uses for DE
    has_raw = full.raw is not None
    log(f"  [{species_id}] cell type '{cell_type}': using {'adata.raw' if has_raw else 'adata.X'} "
        f"as DE matrix", "INFO")
    adata_de = full.raw.to_adata() if has_raw else full.copy()
    adata_de.obs = full.obs.copy()
    adata_de.obs['_fg'] = np.where(fg_mask, 'fg', 'bg')
    adata_de.obs['_fg'] = adata_de.obs['_fg'].astype('category')

    # median sequencing depth of the foreground (audit handle for 'off' calls)
    X = adata_de.X
    fg_counts = np.asarray(X[fg_mask].sum(axis=1)).ravel()
    prov['median_depth'] = float(np.median(fg_counts)) if fg_counts.size else np.nan

    sc.tl.rank_genes_groups(
        adata_de, groupby='_fg', groups=['fg'], reference='bg',
        method=de_method, use_raw=False, pts=True, key_added='de')
    res = sc.get.rank_genes_groups_df(adata_de, group='fg', key='de')

    # the 'off' gate depends on pct_nz_group (from pts=True). If scanpy did not
    # produce it, every frac is NaN and NOTHING will ever be labelled 'off'.
    if 'pct_nz_group' not in res.columns:
        log(f"  [{species_id}] WARNING: 'pct_nz_group' absent from DE result "
            f"(columns: {list(res.columns)}). The 'off' label will never fire for "
            f"this cell type -- check the scanpy version's pts output.", "ERROR")

    # mean expression of foreground per gene, aligned to var_names
    fg_mean = np.asarray(X[fg_mask].mean(axis=0)).ravel()
    mean_by_gene = dict(zip(adata_de.var_names, fg_mean))

    recs = {}
    for _, row in res.iterrows():
        g = row['names']
        frac = row.get('pct_nz_group', np.nan)
        recs[g] = {
            'lfc': float(row['logfoldchanges']),
            'fdr': float(row['pvals_adj']),
            'frac': float(frac) if pd.notna(frac) else np.nan,
            'mean': float(mean_by_gene.get(g, np.nan)),
        }
    log(f"  [{species_id}] cell type '{cell_type}': fg={n_fg} bg={n_bg}, "
        f"{len(recs)} genes tested", "INFO")
    return recs, prov


def label_gene(rec, off_frac, fdr_thr):
    """Assign up / down / flat / off from a DE record (or None if gene absent)."""
    if rec is None:
        return 'off'  # gene not present in this species' matrix at all
    frac = rec['frac']
    if not np.isnan(frac) and frac < off_frac:
        return 'off'
    fdr, lfc = rec['fdr'], rec['lfc']
    if not np.isnan(fdr) and fdr <= fdr_thr and not np.isnan(lfc):
        if lfc > 0:
            return 'up'
        if lfc < 0:
            return 'down'
    return 'flat'


# --------------------------------------------------
def main() -> None:
    args = get_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.input, 'rb') as f:
        sm = pickle.load(f)
    keys = {args.id1: args.anno1, args.id2: args.anno2}
    log(f"Pairwise comparison keys: {keys}", "INFO")

    pms_df = pd.read_csv(args.pms)

    # 1. homology components (BLAST gene graph)
    comps = build_homology_components(sm, args.id1, args.id2, args.blast_thr)
    # id -> gene list, and gene -> component id
    comp_of_gene = {}
    comp_members = {}
    for k, genes in enumerate(comps):
        cid = f"HC_{k:05d}"
        comp_members[cid] = genes
        for g in genes:
            comp_of_gene[g] = cid

    # species of each gene (by prefix)
    species_of_gene = {g: g.split('_', 1)[0] for g in comp_of_gene}

    # Audit: are the component genes actually present in each species' DE matrix?
    # A gene absent here will be silently labelled 'off' (frac 0), which would be
    # a false 'off' driven by matrix pruning rather than biology.
    for sp in (args.id1, args.id2):
        adata = sm.sams[sp].adata
        de_vars = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
        comp_genes_sp = [g for g in comp_of_gene if species_of_gene[g] == sp]
        # strip prefix to compare against var_names, which are unprefixed
        missing = [g for g in comp_genes_sp
                   if (g.split('_', 1)[1] if '_' in g else g) not in de_vars
                   and g not in de_vars]
        log(f"[{sp}] {len(comp_genes_sp)} component genes; {len(missing)} NOT in the DE "
            f"matrix (would be forced to 'off'). Examples: {missing[:10]}", "INFO"
            if not missing else "ERROR")

    # 2. qualifying pairings
    pairings = qualifying_pairings(pms_df, args.id1, args.id2, args.pms_thr)

    # 3. per-(species, cell type) DE, computed once and cached
    side_col = {args.id1: 'celltype_1', args.id2: 'celltype_2'}
    needed = defaultdict(set)  # species -> set(cell_type)
    for pr in pairings:
        needed[args.id1].add(pr['celltype_1'])
        needed[args.id2].add(pr['celltype_2'])

    de_cache = {}     # (species, cell_type) -> gene_records or None
    prov_cache = {}   # (species, cell_type) -> provenance dict
    for sp in (args.id1, args.id2):
        ct_col = keys[sp]
        for ct in sorted(needed[sp]):
            recs, prov = per_celltype_de(sm, sp, ct_col, ct,
                                         args.de_method, args.min_cells)
            de_cache[(sp, ct)] = recs
            prov_cache[(sp, ct)] = prov
            # sanity: how many of THIS species' component genes did we find in
            # the DE result? A near-zero hit rate flags a prefix/name mismatch.
            if recs is not None:
                sp_comp_genes = [g for g in comp_of_gene if species_of_gene[g] == sp]
                hits = sum(1 for g in sp_comp_genes
                           if g in recs or (g.split('_', 1)[1] if '_' in g else g) in recs)
                rate = hits / len(sp_comp_genes) if sp_comp_genes else 0.0
                lvl = "INFO" if rate > 0.5 else "ERROR"
                log(f"  [{sp}] '{ct}': matched {hits}/{len(sp_comp_genes)} component genes "
                    f"to DE result ({rate:.0%})", lvl)

    # 4. walk pairings x components, emit incoherent ones
    rows = []
    for pr in pairings:
        ct1, ct2, pms = pr['celltype_1'], pr['celltype_2'], pr['pms']
        pairing_id = f"{args.id1}:{ct1}__{args.id2}:{ct2}"
        ct_by_species = {args.id1: ct1, args.id2: ct2}

        for cid, genes in comp_members.items():
            member_rows = []
            labels = set()
            for g in genes:
                sp = species_of_gene[g]
                ct = ct_by_species[sp]
                recs = de_cache.get((sp, ct))
                # recs is None => DE was skipped for this cell type (too few
                # cells) => this gene cannot be labelled, so drop it. The
                # component may still be evaluated on its remaining members.
                if recs is None:
                    continue
                # DE 'names' may be species-prefixed ('ax_MMP13') or bare
                # ('MMP13') depending on how var_names survived SAMAP
                # construction. Try both so a prefix mismatch can't silently
                # turn every gene into 'off'.
                bare = g.split('_', 1)[1] if '_' in g else g
                rec = recs.get(g)
                if rec is None:
                    rec = recs.get(bare)
                lab = label_gene(rec, args.off_frac, args.fdr_thr)
                labels.add(lab)
                member_rows.append({
                    'component_id': cid, 'id1': args.id1, 'id2': args.id2,
                    'pairing_id': pairing_id,
                    'celltype_1': ct1, 'celltype_2': ct2, 'pms': pms,
                    'gene': g, 'species': sp, 'celltype': ct, 'label': lab,
                    'log2fc': rec['lfc'] if rec else np.nan,
                    'fdr': rec['fdr'] if rec else np.nan,
                    'frac_expressing': rec['frac'] if rec else 0.0,
                    'mean_expression': rec['mean'] if rec else 0.0,
                })

            # incoherence gate: >1 distinct label among labelled members
            if len(member_rows) >= 2 and len(labels) > 1:
                sig = ','.join(sorted(labels, key=lambda l: ['up', 'down', 'flat', 'off'].index(l)))
                for mr in member_rows:
                    mr['divergence_signature'] = sig
                    mr['n_members'] = len(member_rows)
                rows.extend(member_rows)

    # 5. write outputs
    cols = ['component_id', 'id1', 'id2', 'pairing_id', 'celltype_1', 'celltype_2',
            'pms', 'gene', 'species', 'celltype', 'label', 'log2fc', 'fdr',
            'frac_expressing', 'mean_expression', 'divergence_signature', 'n_members']
    out_main = args.output_dir / f"{args.id1}_{args.id2}_homology_divergence.csv"
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(out_main, index=False)
    if len(df):
        n_comp = df['component_id'].nunique()
        n_pairings_emitted = df['pairing_id'].nunique()
        n_comp_pairing = df.groupby(['component_id', 'pairing_id']).ngroups
        label_counts = df['label'].value_counts().to_dict()
        log(f"Wrote {len(df)} gene rows: {n_comp} distinct components across "
            f"{n_pairings_emitted} pairings, {n_comp_pairing} component x pairing units "
            f"-> {out_main}", "INFO")
        log(f"Label distribution across emitted rows: {label_counts}", "INFO")
        if label_counts.get('off', 0) / len(df) > 0.5:
            log(f"WARNING: >50% of emitted labels are 'off'. If unexpected, check the "
                f"gene-match rates above (possible prefix mismatch) or --off-frac.", "ERROR")
    else:
        log(f"Wrote 0 rows -> {out_main}. No components diverged under any pairing "
            f"(or all were filtered). Check component count and DE match rates above.", "INFO")

    # provenance: one row per (pairing side) actually referenced
    prov_rows = []
    seen = set()
    for pr in pairings:
        for sp, ctkey in ((args.id1, 'celltype_1'), (args.id2, 'celltype_2')):
            ct = pr[ctkey]
            key = (f"{args.id1}:{pr['celltype_1']}__{args.id2}:{pr['celltype_2']}", sp)
            if key in seen:
                continue
            seen.add(key)
            prov = prov_cache.get((sp, ct), {})
            prov_rows.append({
                'pairing_id': key[0], 'species': sp, 'celltype': ct,
                'n_cells': prov.get('n_cells', np.nan),
                'median_depth': prov.get('median_depth', np.nan),
                'n_background_cells': prov.get('n_background_cells', np.nan),
            })
    out_prov = args.output_dir / f"{args.id1}_{args.id2}_pairing_provenance.csv"
    pd.DataFrame(prov_rows).to_csv(out_prov, index=False)
    log(f"Wrote provenance -> {out_prov}", "INFO")


# --------------------------------------------------
if __name__ == '__main__':
    main()
