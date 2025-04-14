// modules/impute_beagle.nf
// Impute genotypes using Beagle 5.5

process BEAGLE5 {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/imputed_beagle5", mode: 'copy'
    

    input:
    tuple val(meta), path(phased_vcf), path(phased_vcf_index), path(ref_panel), path(ref_map_file)
    val(ch_java_memory) // Memory for java task

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}.vcf.gz"), emit: imputed_vcf


    script:
    prefix = "${meta.id}_imputed_${meta.target_cohort}"
    
    """
    echo "Starting imputation for ${meta.id} (Target: ${meta.target_cohort}) using ${task.cpus} CPUs"
    echo "Reference panel file: ${ref_panel}"
    echo "Input VCF: ${phased_vcf}"

    java -jar -Xmx${ch_java_memory}G  \
        ${params.beagle5} \
        gt=${phased_vcf} \
        map=${ref_map_file} \
        ref=${ref_panel} \
        out=${prefix} \
        gp=true \
        ap=true \
        nthreads=${task.cpus}
        
    echo "Imputation completed for ${meta.id} with reference panel from ${meta.target_cohort}"
    """
}