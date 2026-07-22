/*
 *  MODULE: homology_component_de.nf
 *
 *  Description:
 *      Per-cell-type divergence of BLAST-linked (homology) gene components.
 *      For each cross-species cell-type pairing with PMS >= threshold, walks
 *      each connected component of the BLAST gene graph (sm.gnnm) and returns
 *      the component iff its member genes do not all share one expression
 *      label {up, down, flat, off} on their own species' side of the pairing.
 *
 *      NOTE: a "homology component" here is a connected component of the BLAST
 *      *gene* graph -- this is NOT the cell-type "alignment family" concept
 *      used in connected_de.nf. Different graph, different object.
 *
 *      This module is a terminal sink: it depends only on SUMMARY_SAMAP
 *      outputs and nothing downstream consumes its output. It runs in
 *      parallel with CONNECTED_DE.
 *
 *  Inputs:
 *      run_id:        Run identifier
 *      tuple(id1, id2, anno1, anno2, pms, samap_obj):
 *          id1/id2:       the pairwise species-species comparison
 *          anno1/anno2:   annotation (obs) columns for id1/id2
 *          pms:           pms_cluster_alignment_scores.csv from SUMMARY_SAMAP
 *          samap_obj:     cleaned SAMap object (samap_results_Cleaned.pkl)
 *
 *  Outputs:
 *      <id1>_<id2>_homology_divergence.csv    tidy, one row per (component x pairing x gene)
 *      <id1>_<id2>_pairing_provenance.csv     one row per (pairing x species side)
 *      logfile
 */

process HOMOLOGY_COMPONENT_DE {
    tag "${run_id} - Per-cell-type Homology Component Divergence"

    container 'mdiblbiocore/postanalysis:latest'

    input:
        val run_id
        tuple val(id1), val(id2), val(anno1), val(anno2), path(pms), path(samap_obj)

    output:
        tuple val(id1), val(id2), path("${id1}_${id2}_homology_divergence.csv"), emit: divergence
        tuple val(id1), val(id2), path("${id1}_${id2}_pairing_provenance.csv"),  emit: provenance
        path "${run_id}_homologyComponentDE.log"

    script:
    """
    LOG="${run_id}_homologyComponentDE.log"
    homology_component_de.py \\
        --input ${samap_obj} \\
        --pms ${pms} \\
        --id1 ${id1} --anno1 ${anno1} \\
        --id2 ${id2} --anno2 ${anno2} \\
        --output_dir . 2>&1 | tee -a \$LOG
    """
}
