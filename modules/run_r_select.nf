#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// modules/run_r_select.nf
// Runs the SNP_Selection.R script to generate filter lists

process RUN_R_SELECT {
    tag "$imputed_cohort"
    publishDir "${params.outdir}/${imputed_cohort}/R_output", mode: 'copy'

    
    input:
    tuple val(imputed_cohort), val(target_cohort_gather), path(info_summary_file) // e.g., ["cohortA", "path/to/cohortA_INFO_files.txt"]

    output:
    // Use file glob patterns that match potential outputs of the R script
    tuple val(imputed_cohort), val(target_cohort_gather), path("Selected_filterout*.txt"), emit: filterlists
    tuple val(imputed_cohort), val(target_cohort_gather), path("Selected_keeplist*.txt"), optional: true, emit: keeplists
    tuple val(imputed_cohort), val(target_cohort_gather), path("Summary_keeplist*.txt"), optional: true, emit: sum_keeplist
    tuple val(imputed_cohort), val(target_cohort_gather), path("Summary_filterout*.txt"), optional: true, emit: sum_filterout

    script:
    // The R script is expected to run in the current directory and take the info summary file as input
    // It should generate output files (keeplist_*.txt or filterout_*.txt) in the current directory.
    """
    #!/bin/bash
    set -euo pipefail

    echo "Running R script ${params.snp_select_r} on ${info_summary_file}"
    Rscript ${params.snp_select_r} ${info_summary_file}

    # Check if at least one filter or keep list was generated
    # List files matching the patterns and check if the list is non-empty
    filter_count=\$(ls Selected_filterout*.txt 2>/dev/null | wc -l || echo 0)
    keep_count=\$(ls Selected_keeplist*.txt 2>/dev/null | wc -l || echo 0)
    
    # Debug: List all files in the current directory
    echo "Files generated in the current directory:"
    ls -l
    

    if [ \$filter_count -eq 0 ] && [ \$keep_count -eq 0 ]; then
      echo "Warning: R script '${params.snp_select_r}' did not generate any Selected_filterout*.txt or Selected_keeplist*.txt files for cohort ${imputed_cohort}." >&2
      echo "Creating dummy files to prevent pipeline failure." >&2
      touch Selected_filterout_dummy_${imputed_cohort}.txt
      touch Selected_keeplist_dummy_${imputed_cohort}.txt
    else
       echo "R script generated \$filter_count filterout file(s) and \$keep_count keeplist file(s)."
    fi
    """
}
