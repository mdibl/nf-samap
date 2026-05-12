include { CREATE_LOUPE_INPUT } from '../modules/create_loupe_input.nf'
include { CREATE_LOUPE_FILE } from '../modules/create_loupe_file.nf'


workflow CREATE_LOUPE {
    take:
        run_id
        compTuple       // channel tuple containing everything necessary to build Loupe Object

    main:
        CREATE_LOUPE_INPUT(
                run_id,
                compTuple.map { id1, id2, anno1, anno2, h5ad1, h5ad2, _samap -> [id1, id2, anno1, anno2, h5ad1, h5ad2] },
                compTuple.map { _id1, _id2, _anno1, _anno2, _h5ad1, _h5ad2, samap -> samap }.first()
        )
        CREATE_LOUPE_FILE(
            run_id,
            CREATE_LOUPE_INPUT.out.loupeinput
        )

}