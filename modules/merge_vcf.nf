// Process 9: MERGE_VCF
// Merge all filtered VCF files for each chromosome
process MERGE_VCF {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}/merged", mode: 'copy'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}_combined.vcf.gz"), path("${meta.id}_combined.vcf.gz.tbi"), emit: merged_vcf

    script:
    def output_prefix = "${meta.id}_combined"
    // Create a space-separated list of VCF files
    def vcf_list = vcfs.join(" ")
    """
    echo "Merging filtered VCF files for ${meta.id}"
    
    # Merge VCF files using bcftools concat
    ${params.bcftools} concat -a ${vcf_list} -Oz -o ${output_prefix}.vcf.gz
    
    # Index the merged output
    ${params.tabix} -p vcf ${output_prefix}.vcf.gz
    
    echo "VCF merging completed for ${meta.id}"
    """
} 