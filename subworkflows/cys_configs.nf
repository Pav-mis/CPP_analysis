include { GET_CYS_CONFIGS; CLUSTER_CONFIGS} from "../modules/modules.nf"

workflow CYS_CONFIGS {
    take:
    intermediate

    main:
    GET_CYS_CONFIGS(intermediate)
    CLUSTER_CONFIGS(GET_CYS_CONFIGS.out)

    emit:
    configs = GET_CYS_CONFIGS.out
    positions = CLUSTER_CONFIGS.out
}