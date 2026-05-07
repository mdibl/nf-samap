include { CREATE_LOUPE_INPUT } from '../modules/create_loupe_input.nf'
include { CREATE_LOUPE_FILE } from '../modules/create_loupe_file.nf'


workflow CREATE_LOUPE {
    take:
        run_id
        compTuple       // channel tuple containing everything necessary to build Loupe Object

    main:
        CREATE_LOUPE_INPUT(
            run_id,
            compTuple
        )
        CREATE_LOUPE_FILE(
            run_id,
            CREATE_LOUPE_INPUT.out.loupeinput
        )

}