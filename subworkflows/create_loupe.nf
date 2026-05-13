include { CREATE_LOUPE_INPUT } from '../modules/create_loupe_input.nf'
include { CREATE_LOUPE } from       '../modules/create_loupe.nf'

workflow CREATE_LOUPE_WORKFLOW {
    take:
        run_id
        compTuple   // [id1, id2, anno1, anno2, h5ad1, h5ad2, samap_cleaned, alignment_families]

    main:
        CREATE_LOUPE_INPUT(
            run_id,
            compTuple.map { id1, id2, anno1, anno2, h5ad1, h5ad2, _samap, _af -> 
                [id1, id2, anno1, anno2, h5ad1, h5ad2] 
            },
            compTuple.map { _id1, _id2, _anno1, _anno2, _h5ad1, _h5ad2, samap, _af -> 
                samap 
            }.first(),
            compTuple.map { _id1, _id2, _anno1, _anno2, _h5ad1, _h5ad2, _samap, af -> 
                af 
            }
        )

        CREATE_LOUPE(
            run_id,
            CREATE_LOUPE_INPUT.out.loupeinput
        )
}