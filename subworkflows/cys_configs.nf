include { GET_CYS_CONFIGS; CLUSTER_CONFIGS} from "../modules/modules.nf"

workflow CYS_CONFIGS {
    take:
    intermediate

    main:
    GET_CYS_CONFIGS(intermediate)
    CLUSTER_CONFIGS(GET_CYS_CONFIGS.out.configs)

    emit:
    configs = GET_CYS_CONFIGS.out.configs
    positions = CLUSTER_CONFIGS.out
    fasta = GET_CYS_CONFIGS.out.fasta
}