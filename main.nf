include { CYS_CONFIGS} from "./subworkflows/cys_configs.nf"
include { GET_CPP_SITES} from "./subworkflows/get_cpp_sites.nf"
include { GEN_INTERMEDIATE} from "./subworkflows/generate_intermediate.nf"
include { PLOT_CONFIGS; BED_TO_TABLE; CALC_HYDROPHOBICITY; BED_TO_TABLE_WIDE; CREATE_SEQ_TABLES; PLOT_SEQ_TABLE; GET_DREK_MERS} from "./modules/modules.nf"

workflow {
    protiens = file(params.protiens)
    cpp_file = file(params.cpp_file)
    pred_ranked = file(params.pred_ranked)
    pred_gff = file(params.pred_gff)

    
    GEN_INTERMEDIATE(protiens, pred_ranked, pred_gff)
    CYS_CONFIGS(GEN_INTERMEDIATE.out)
    configs = CYS_CONFIGS.out.positions.flatten()
    CALC_HYDROPHOBICITY(protiens)
    GET_DREK_MERS(protiens)
    GET_CPP_SITES(cpp_file, protiens)
    BED_TO_TABLE(GET_CPP_SITES.out, CYS_CONFIGS.out.configs)
    BED_TO_TABLE_WIDE(GET_CPP_SITES.out, GEN_INTERMEDIATE.out, GET_DREK_MERS.out)
    tables = CREATE_SEQ_TABLES(BED_TO_TABLE_WIDE.out, CALC_HYDROPHOBICITY.out).flatten()
    plots = PLOT_SEQ_TABLE(tables)
    
    //plots = PLOT_CONFIGS(configs, BED_TO_TABLE.out, GEN_INTERMEDIATE.out)
}