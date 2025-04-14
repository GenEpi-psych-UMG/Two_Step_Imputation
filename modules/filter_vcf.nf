#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// modules/filter_vcf.nf
// Filters imputed VCFs based on keeplist or filterout list

process FILTER_VCF {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/filtered_vcfs", mode: 'copy'  
                    
    input:
    tuple val(meta), path(imputed_vcf), path(filter_list) // meta e.g., [id:"cohortA_chr1_imputed_cohortB", cohort:"cohortA", ref:"cohortB", chr:1]
    // filter_list is EITHER keeplist_cohortB.txt OR filterout_cohortC.txt (or similar)
    
    output:    
    tuple val(meta), path("*.recode.vcf.gz"), emit: filtered_vcf    
    tuple val(meta), path("*.recode.vcf.gz.tbi"), optional:true, emit: filtered_idx    
        
    script:
    // Define filter logic here where task context is available
    def filter_list_name = filter_list.baseName
    def filter_option = filter_list_name.startsWith('Selected_keeplist') ? "--positions ${filter_list}" : "--exclude-positions ${filter_list}"
    def filter_type_suffix = filter_list_name.startsWith('Selected_keeplist') ? 'addon' : 'kept'
    def out_filename = "${meta.id}_filtered_${meta.target_cohort}_${filter_type_suffix}.recode.vcf.gz"
    
    """
    #!/bin/bash
    set -euo pipefail
    echo "Filtering ${imputed_vcf} using VCFtools with option: ${filter_option}"
    
    if [ ! -s "${filter_list}" ]; then        
        echo "Warning: Filter list '${filter_list}' is empty or not found for ${meta.id}."    
        
    else
        vcftools \
        --gzvcf ${imputed_vcf} \
        ${filter_option} \
        --recode \
        --recode-INFO-all \
        --stdout | ${params.bgzip} -c > ${out_filename}
        
        tabix -p vcf ${out_filename}    
    fi    
    """
}
