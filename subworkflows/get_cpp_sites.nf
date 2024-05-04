include { BLAST_DB; BLASTP; BLAST_TO_BED; SORT_AND_MERGE_BED} from "../modules/modules.nf"

workflow GET_CPP_SITES {
    take:
    cpp_file
    protiens

    main:
    BLAST_DB(protiens)
    BLASTP(protiens, BLAST_DB.out, cpp_file)
    BLAST_TO_BED(BLASTP.out)
    SORT_AND_MERGE_BED(BLAST_TO_BED.out)

    emit:
    SORT_AND_MERGE_BED.out
}