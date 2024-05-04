include { FASTA_TO_TABLE; GENERATE_INTERMEDIATE} from "../modules/modules.nf"

workflow GEN_INTERMEDIATE {
    take:
    protiens
    pred_ranked
    pred_gff

    main:
    FASTA_TO_TABLE(protiens)
    GENERATE_INTERMEDIATE(pred_ranked, pred_gff, FASTA_TO_TABLE.out)

    emit:
    GENERATE_INTERMEDIATE.out
}