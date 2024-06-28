process FASTA_TO_TABLE {
    container 'docker://quay.io/biocontainers/biopython:1.81'
    
    input:
    path protiens

    output:
    path "${protiens.baseName}.tsv"

    """
    python3 $projectDir/bin/fasta2table.py ${protiens} ${protiens.baseName}.tsv
    """
}

process GENERATE_INTERMEDIATE {
    container 'docker://amancevice/pandas:2.2.2'
    publishDir "${params.outDir}"

    input:
    path pred_ranked
    path pred_gff
    path protien_table

    output:
    path "intermediate.tsv"

    """
    python3 $projectDir/bin/createIntermediate.py ${pred_ranked} ${pred_gff} ${protien_table}
    """
}

process GET_CYS_CONFIGS {
    input: 
    path protiens

    output:
    path "${protiens.baseName}.cys.txt", emit: configs
    path "${protiens.baseName}_mature.fasta", emit: fasta

    """
    perl $projectDir/bin/cys_classify.pl ${protiens} ${protiens.baseName}.cys.txt
    perl $projectDir/bin/mature_seqs.pl ${protiens} ${protiens.baseName}_mature.fasta
    """
}

process CLUSTER_CONFIGS {
    container 'docker://quay.io/biocontainers/pandas:1.5.2'

    input:
    path cys

    output:
    path("*.positions.txt")


    """
    python3 $projectDir/bin/clusterConfigs.py ${cys}
    """
}

process PLOT_CONFIGS {
    container 'docker://quay.io/uphl/seaborn:0.12.2-2-slim'
    publishDir "${params.outDir}/cys_configs"
    cache false


    input:
    path positions
    path cpp
    path intermediate


    output:
    path "${positions.baseName}.png"
    path "${positions.baseName}.intermediate.tsv"

    """
    python3 $projectDir/bin/plotConfigs.py '${positions}' '${cpp}' '${intermediate}' '${positions.baseName}'
    """


}

process BLAST_DB {
    container 'docker://quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    input:
    path protiens

    output:
    path "${protiens.baseName}*"

    """
    makeblastdb -in ${protiens} -dbtype prot
    """

}

process BLASTP {
    container 'docker://quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    input:
    path protiens
    path db
    path cpp_file

    output:
    path "CPP_sites.out"

    """
    blastp -db ${protiens} -query ${cpp_file} -evalue 1000 -word_size 2 -matrix PAM30 -gapopen 10 -gapextend 1 -seg no -out CPP_sites.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs"
    """

}

process BLAST_TO_BED {
    container 'docker://quay.io/biocontainers/pandas:1.5.2'

    input:
    path blast

    output:
    path "CPP_sites.bed"

    """
    python3 $projectDir/bin/blast2bed.py ${blast}
    """


}



process SORT_AND_MERGE_BED {
    container "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"

    input:
    path bed

    output:
    path "CPP_sites.sorted.merged.bed"

    """
    bedtools sort -i ${bed} > CPP_sites.sorted.bed
    bedtools merge -i CPP_sites.sorted.bed > CPP_sites.sorted.merged.bed
    """
}

process BED_TO_TABLE{
    container 'docker://quay.io/biocontainers/pandas:1.5.2'

    input:
    path bed
    path cys

    output:
    path "CPP_sites.positions.txt"

    """
    python3 $projectDir/bin/bed2table.py ${bed} ${cys}
    """
}

process BED_TO_TABLE_WIDE{

    container 'docker://quay.io/biocontainers/pandas:1.5.2'

    input:
    path bed
    path cys
    path drek

    output:
    path "CPP_sites.positions.txt"

    """
    python3 $projectDir/bin/bed2table_wide.py ${bed} ${cys} ${drek}
    """
}

process CALC_HYDROPHOBICITY{
    container 'docker://quay.io/biocontainers/emboss:6.6.0--hdde3b0b_8'
    publishDir "${params.outDir}"

    input:
    path mature_seqs

    output:
    path "hydropathy/*"

    """
    mkdir -p split_sequences
    mkdir -p hydropathy

    awk 'BEGIN {output_file = ""} 
         /^>/ {
             if (output_file != "") close(output_file)
             header = substr(\$0, 2)
             gsub(/ /, "_", header)
             gsub(/[^a-zA-Z0-9_#.,()+-]/, "", header)
             output_file = "split_sequences/" header ".fasta"
             print \$0 > output_file
         }
         /^[^>]/ { print \$0 >> output_file }
        ' ${mature_seqs}

    for file in split_sequences/*.fasta; do
        line_count=\$(wc -l < "\$file")
        if [ "\$line_count" -gt 1 ]; then
            second_line=\$(sed -n '2p' "\$file")
            second_line_length=\$(echo -n "\$second_line" | wc -m)
            if [ "\$second_line_length" -gt 1 ]; then
                base_name=\$(basename "\$file" .fasta)
                pepwindowall -sequence "\$file" -graph data -window 21 -goutfile "\$base_name"
            fi
        fi
    done

    mv *.dat hydropathy
    
    for file in "hydropathy"/*
    do
        if [ -f "\$file" ]; then
    
            grep -v '^#' "\$file" > temp_file && mv temp_file "\$file"
        fi
    done
    
    """


}

process CREATE_SEQ_TABLES {
    container 'docker://quay.io/biocontainers/pandas:1.5.2'

    input:
    path table_wide
    path hydropathy

    output:
    path "tables/*"

    """
    mkdir hydropathy
    mv *.dat hydropathy

    mkdir tables
    python3 $projectDir/bin/seqtable.py ${table_wide} hydropathy tables
    """



}

process PLOT_SEQ_TABLE {

    cache false
    publishDir "/scratch/y95/pmisiun/hydro_plots", mode: 'copy'
    container 'docker://joseespinosa/docker-r-ggplot2'

    input: 
    path seq_table

    output:
    path "${seq_table.baseName}.png"

    """
    Rscript $projectDir/bin/plot_curve.r  '${seq_table}' '${seq_table.baseName}'
    """
}

process GET_DREK_MERS {

    input:
    path fasta

    output:
    path "${fasta.baseName}.drek.tsv"

    """
    $projectDir/bin/DREK_mers.pl ${fasta} ${fasta.baseName}.drek.tsv
    """


}

