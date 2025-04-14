#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// modules/extract_info.nf
// Runs the extract_info_VCF.pl script to generate .all.info.txt files

process EXTRACT_INFO {
    tag "($meta.id (target:$meta.target_cohort))" 
    publishDir "${params.outdir}/${meta.cohort}/info_scores", mode: 'copy', pattern: '*.all.info.txt'
  
    input:
    tuple val(meta), path(imputed_vcf) 

    output:
    tuple val(meta), path("*.all.info.txt"), emit: info_file

    // Note: \$ refers to bash variables, normal $ is for nextflow variables
    script:
    
    // Choose extract info script based on imputer
    def extract_info_script = params.imputer == "beagle5" ? params.extract_info_pl_beagle5 : params.extract_info_pl

    """
    #!/bin/bash
    set -euo pipefail

    baseName=\$(basename ${imputed_vcf} .vcf.gz)
  
    perl ${extract_info_script} ./ \$baseName

    # Check if the expected output file exists
    if [ ! -f "\${baseName}.all.info.txt" ]; then
      echo "Error: Perl script '${extract_info_script}' did not generate expected output \${baseName}.all.info.txt" >&2
      exit 1
    fi
    """
}