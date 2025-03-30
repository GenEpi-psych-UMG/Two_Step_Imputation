// Process 7: RUN_R_SELECT
// Run R script to analyze info files and generate filter lists
process RUN_R_SELECT {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}/r_output", mode: 'copy'
    
    conda "r-base=4.1.0"

    input:
    tuple val(meta), path(info_files)

    output:
    tuple val(meta), path("Keep_list_main_chr*.txt"), emit: filter_lists
    
    script:
    def info_files_list = "info_files_list.txt"
    """
    echo "Running R SNP selection script for ${meta.cohort}"
    
    # Create the input file listing all info files
    touch ${info_files_list}
    for file in ${info_files}; do
        echo "\${file}" >> ${info_files_list}
    done
    
    # Run the R script
    Rscript ${params.snp_select_r} ${info_files_list}
    
    echo "R SNP selection completed for ${meta.cohort}"
    """
} 