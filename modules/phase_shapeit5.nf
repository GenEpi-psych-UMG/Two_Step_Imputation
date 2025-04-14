// modules/phase_shapeit5.nf
// Phase genotypes using SHAPEIT5 (Phase stage)

process SHAPEIT5 {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}/phased_shapeit5", mode: 'copy'

    // SHAPEIT5 needs bcftools for input conversion/prep if needed
    // Assuming SHAPEIT5 path is handled by beforeScript in config

    input:
    tuple val(meta), path(bcf_file), path(bcf_index), path(ref_map_chr)

    output:
    tuple val(meta), path("${meta.id}_phased_shapeit5.bcf"), path("${meta.id}_phased_shapeit5.bcf.csi"), emit: phased_bcf // SHAPEIT5 outputs BCF

    script:
    // Define output_prefix for use within the script commands
    def output_prefix = "${meta.id}_phased_shapeit5"

    """
    echo "Starting phasing for ${meta.id} with SHAPEIT5 using ${task.cpus} threads"

    ${params.bcftools} --version

    # SHAPEIT5 Phase command
    ${params.shapeit5} \
        --input ${bcf_file} \
        --region ${meta.chr} \
        --map ${ref_map_chr} \
        --output ${output_prefix}.bcf \
        --output-format bcf \
        --thread ${task.cpus} \
        --log ${output_prefix}.log 

     ${params.bcftools} index -c -f ${output_prefix}.bcf


    echo "SHAPEIT5 phasing completed for ${meta.id}"
    """
} 