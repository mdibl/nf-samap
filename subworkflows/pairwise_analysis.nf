include { SUMMARY_SAMAP } from '../modules/summary_samap.nf'
include { CONNECTED_DE } from '../modules/connected_de.nf'
include { ADDITIONAL_ANALYSIS } from '../modules/additional_analysis.nf'

workflow PAIRWISE_ANALYSIS {
    take:
        run_id       // val run_id
        ch_samples      // [meta, h5ad, fasta]
        samap_results   // RUN_SAMAP.out.results


    main:
        // Necessary Context
        id_anno = ch_samples
            .map { meta, _h5ad, _fasta -> [meta.id, meta.annotation] }

        //Make id + Annotation pairs for pairwise post-analysis
        idCompare = id_anno.combine(id_anno)
            .filter { a_id, _a_anno, b_id, _b_anno -> a_id < b_id }
            .map { a_id, a_anno, b_id, b_anno ->
                [a_id, b_id, a_anno, b_anno]
        }

        id_anno_collected = id_anno
            .collect()
            .map { pairs -> 
                [pairs.collate(2).collect { p -> p[0] }, pairs.collate(2).collect { r -> r[1] }]
            }


        SUMMARY_SAMAP(
            run_id,
            idCompare.combine(samap_results),
            id_anno_collected
        )

        CONNECTED_DE(
            run_id,
            idCompare
                .join(SUMMARY_SAMAP.out.pms,          by: [0, 1]) // match on id1, id2
                .join(SUMMARY_SAMAP.out.samap_cleaned, by: [0, 1]) // match on id1, id2
        )


        // To line up corect order of pairwise comparison with correct upstream output
        combined = CONNECTED_DE.out.groupinganalysis
            .join(CONNECTED_DE.out.analysis,        by: [0,1])
            .join(CONNECTED_DE.out.all_de_results,  by: [0,1])
            .join(SUMMARY_SAMAP.out.genepairs,      by: [0,1])
            .join(SUMMARY_SAMAP.out.pms,            by: [0,1])
            .join(idCompare,                        by: [0,1])


        ADDITIONAL_ANALYSIS(
            run_id,
            combined
        )

    //emit:
        //results = ADDITIONAL_ANALYSIS.out.results
}