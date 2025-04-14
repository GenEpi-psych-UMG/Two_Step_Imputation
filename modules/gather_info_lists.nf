#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// modules/gather_info_lists.nf
// Gathers all .all.info.txt file paths for a single *imputed* cohort
// and creates a summary file (e.g., CohortA_INFO_files.txt)

process GATHER_INFO_LISTS {
    tag "$cohort"
    publishDir "${params.outdir}/${cohort}/R_input", mode: 'copy'

    input:
    tuple val(cohort), val(target_cohort_gather), path(info_files) // cohort is the cohort being imputed (e.g., "cohortA")
    val ordered_cohorts

    output:
    tuple val(cohort), val(target_cohort_gather), path("${cohort}_INFO_files.txt"), emit: info_summary

    script:
    def info_list_filename = "${cohort}_INFO_files.txt"
    def temp_file = "${info_list_filename}.temp"
    """
    #!/bin/bash
    set -euo pipefail

    echo -e "Path\tReference\tchr" > ${temp_file}

    # Create a temporary file to store the entries
    for info_file in ${info_files.join(' ')}; do
        # Extract Reference cohort and Chromosome from the info file path/name
        # Filename format like: /path/to/imputedCohort_chr96_imputed_refCohort.all.info.txt
        base=\$(basename "\$info_file" .all.info.txt)
        ref_cohort=\$(echo "\$base" | sed -n 's/.*_imputed_\\(.*\\)/\\1/p')
        chr=\$(echo "\$base" | sed -n 's/.*_chr\\([0-9]*\\)_imputed_.*/\\1/p')

        # Check if extraction was successful
        if [ -z "\$ref_cohort" ] || [ -z "\$chr" ]; then
            echo "Warning: Could not extract reference cohort or chromosome from \$info_file. Skipping." >&2
            continue
        fi

        # Get the absolute path to the info file
        absolute_info_file_path=\$(readlink -f \$info_file)

        echo -e "\$absolute_info_file_path\t\$ref_cohort\t\$chr" >> ${temp_file}
    done

    # Create a temporary file with cohort order
    echo "${ordered_cohorts.join('\n')}" > cohort_order.txt

    # Sort the file by reference cohort (using the order from cohort_order.txt) and then by chromosome
    # Script get associated order 1 2 3 from cohort_order.txt and then appends it to the temp file
    # The number on the 4th column is organized by the order of the 2nd column (reference cohort)
    # Additionally, the 3rd column (chromosome) is sorted numerically
    # This step ensures that the SNP_Selection stepp will process in the correct order!!!!!
    awk -F'\t' '
    NR==FNR { order[\$0]=NR; next }
    FNR==1 { print; next }
    { print \$0 "\t" order[\$2] }
    ' cohort_order.txt ${temp_file} | sort -k4,4n -k3,3n | \
    cut -f1-3 > ${info_list_filename}
    
    # Clean up temporary files
    rm ${temp_file} cohort_order.txt
    """
}