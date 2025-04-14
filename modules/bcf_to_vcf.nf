// modules/bcf_to_vcf.nf
process BCF_TO_VCF {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}/phased_vcf_converted", mode: 'copy'

    // REMOVED Conda directive - Relying on params.bcftools and params.tabix paths
    // Ensure these tools are accessible via PATH (e.g., set in config's process block)
    // Example config addition:
    // withName: 'BCF_TO_VCF' {
    //     beforeScript = 'export PATH="${file(params.bcftools).parent}:${file(params.tabix).parent}:$PATH"'
    // }

    input:
    tuple val(meta), path(bcf_file), path(bcf_index)

    output:
    tuple val(meta), path("${meta.id}_phased_shapeit5.vcf.gz"), path("${meta.id}_phased_shapeit5.vcf.gz.tbi"), emit: converted_vcf

    script:
    def output_vcf = "${meta.id}_phased_shapeit5.vcf.gz"
    """
    echo "Converting BCF ${bcf_file} to VCF.gz for ${meta.id}"
    bcftools --version
    tabix --version || true

    bcftools view \
        --threads ${task.cpus} \
        -O z \
        -o ${output_vcf} \
        ${bcf_file}

    if [ -s "${output_vcf}" ]; then
        tabix -p vcf ${output_vcf}
    else
        echo "ERROR: BCF to VCF conversion failed or produced empty file: ${output_vcf}" >&2
        exit 1
    fi
    echo "BCF to VCF conversion complete for ${meta.id}"
    """
} 