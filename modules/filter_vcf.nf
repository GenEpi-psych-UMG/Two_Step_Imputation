// Process 8: FILTER_VCF
// Filter imputed VCF files based on R output
process FILTER_VCF {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/filtered/${meta.target_cohort}", mode: 'copy'
    
    conda "bioconda::vcftools=0.1.17"

    input:
    tuple val(meta), path(imputed_vcf), path(filter_lists)

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}_filtered.vcf.gz"), emit: filtered_vcf

    script:
    def output_prefix = "${meta.id}_imputed_${meta.target_cohort}_filtered"
    def chr = meta.chr
    def filter_file = "Keep_list_main_chr${chr}_${meta.target_cohort}.txt"
    """
    echo "Filtering ${imputed_vcf} using ${filter_file}"
    
    # Filter VCF using the filter list from R
    if [ -f "${filter_file}" ]; then
        vcftools --gzvcf ${imputed_vcf} \
                 --chr ${chr} \
                 --positions ${filter_file} \
                 --recode \
                 --out ${output_prefix}
    else
        # Alternative filter if specific list doesn't exist
        vcftools --gzvcf ${imputed_vcf} \
                 --chr ${chr} \
                 --recode \
                 --out ${output_prefix}
        echo "Warning: ${filter_file} not found, using all variants"
    fi
    
    # Compress and index the output
    bgzip -c ${output_prefix}.recode.vcf > ${output_prefix}.vcf.gz
    ${params.tabix} -p vcf ${output_prefix}.vcf.gz
    
    echo "VCF filtering completed for ${meta.id}"
    """
} 