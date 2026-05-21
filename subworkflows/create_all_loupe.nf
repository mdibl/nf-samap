include { CREATE_ALL_LOUPE_INPUT } from '../modules/create_all_loupe_input.nf'
include { CREATE_ALL_LOUPE_FILE  } from '../modules/create_all_loupe_file.nf'

workflow CREATE_ALL_LOUPE_WORKFLOW {
    take:
        run_id
        all_info   // channel emitting [ids, annos, h5ads] — collected lists for all species
        samap      // SAMap pickle (RUN_SAMAP.out.results, uncleaned)

    main:
        CREATE_ALL_LOUPE_INPUT(
            run_id,
            all_info,
            samap
        )

        CREATE_ALL_LOUPE_FILE(
            run_id,
            CREATE_ALL_LOUPE_INPUT.out.loupeinput
        )
}
