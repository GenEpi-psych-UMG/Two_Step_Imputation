// Process 5: IMPUTE
// Impute genotypes using minimac4
process IMPUTE {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/imputed/${meta.target_cohort}", mode: 'copy'
    
    conda "bioconda::minimac4=4.1.4"

    input:
    tuple val(meta), path(phased_vcf), path(phased_vcf_index)
    val ref_panels_map

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}.vcf.gz"), emit: imputed_vcf

    script:
    def target_key = [chr: meta.chr, cohort: meta.target_cohort]
    def ref_panel = ref_panels_map.get(target_key)
    def output_prefix = "${meta.id}_imputed_${meta.target_cohort}"
    
    if (!ref_panel) {
        error "Reference panel not found for chromosome ${meta.chr} and cohort ${meta.target_cohort}"
    }
    
    """
    echo "Imputing ${meta.id} using reference panel from ${meta.target_cohort}"
    
    # Run minimac4 imputation
    minimac4 \
      ${ref_panel} \
      ${phased_vcf} \
      --output ${output_prefix}.vcf.gz \
      -f GT \
      --all-typed-sites \
      --threads ${task.cpus}
    
    echo "Imputation completed for ${meta.id} with reference panel from ${meta.target_cohort}"
    """
} 