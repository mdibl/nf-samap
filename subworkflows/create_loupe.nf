include { CREATE_LOUPE_INPUT } from '../modules/create_loupe_input.nf'
include { CONNECTED_DE } from '../modules/connected_de.nf'


workflow CREATE_LOUPE {
    take:
        run_id
        compTuple       // channel tuple containing everything necessary to build Loupe Object
        
    main:
        CREATE_LOUPE_INPUT(
            run_id,
            compTuple
        )

    //emit:
       // pairCompare = idCompare
}