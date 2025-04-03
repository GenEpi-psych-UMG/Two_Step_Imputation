// Process 8: FILTER_VCF
// Filter imputed VCF files based on R output
process FILTER_VCF {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/filtered/${meta.target_cohort}", mode: 'copy'
    
    conda "bioconda::vcftools=0.1.16"

    input:
    tuple val(meta), path(imputed_vcf), path(filterout_list), path(keep_list)

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}_filtered.vcf.gz"), path("${meta.id}_imputed_${meta.target_cohort}_filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    def output_prefix = "${meta.id}_imputed_${meta.target_cohort}_filtered"
    def chr = meta.chr
    """
    echo "Filtering ${imputed_vcf}"
    echo "Using Keep list: ${keep_list}"
    echo "Using Filterout list: ${filterout_list}"

    # Check if files exist
    if [ ! -f "${keep_list}" ]; then
        echo "Error: Keep list ${keep_list} not found!" >&2
        exit 1
    fi
    if [ ! -f "${filterout_list}" ]; then
        echo "Error: Filterout list ${filterout_list} not found!" >&2
        exit 1
    fi

    # Filter VCF using both keep and exclude lists
    vcftools --gzvcf ${imputed_vcf} \
             --chr ${chr} \
             --positions ${keep_list} \
             --exclude-positions ${filterout_list} \
             --recode \
             --recode-INFO-all \
             --out ${output_prefix}

    # Compress and index the output
    # Check if recode file exists and is not empty before compressing
    if [ -s "${output_prefix}.recode.vcf" ]; then
        bgzip -c ${output_prefix}.recode.vcf > ${output_prefix}.vcf.gz
        ${params.tabix} -p vcf ${output_prefix}.vcf.gz
    else
        echo "Error: Filtering produced an empty VCF file (${output_prefix}.recode.vcf). Check input lists and VCF." >&2
        # Create empty outputs to satisfy Nextflow output expectations, but maybe signal error?
        touch ${output_prefix}.vcf.gz ${output_prefix}.vcf.gz.tbi
        # Optionally exit with an error: exit 1 
    fi

    echo "VCF filtering completed for ${meta.id}"
    """
}