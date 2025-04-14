#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// modules/merge_vcf.nf
// Merges filtered VCF files for a given cohort and chromosome

process MERGE_VCF {
    // Use task.ext to make cohort/chr available outside script block
    // beforeScript '''
    // task.ext.cohort = cohort_chr_id[0]
    // task.ext.chr = cohort_chr_id[1]
    // '''

    tag "${meta.id}(${filtered_vcfs.size()} paritions)"
    publishDir "${params.outdir}/${meta.cohort}/merged_vcfs", mode: 'copy'

    input:
    tuple val(meta), path(filtered_vcfs)

    output:
    tuple val(meta), path("${meta.id}.merged.vcf.gz"), emit: merged_vcf
    tuple val(meta), path("${meta.id}.merged.vcf.gz.tbi"), emit: merged_idx

    script:
    // Use task.ext here as well
    def file_list = "files_to_merge.txt"
    def out_vcf = "${meta.id}.merged.vcf.gz"

    """
    #!/bin/bash
    set -euo pipefail

    echo "Merging ${filtered_vcfs.size()} VCF(s) for ${meta.cohort} chromosome ${meta.chr}"

    # Create a file listing the VCFs to merge
    for vcf_file in ${filtered_vcfs.join(' ')}; do
        readlink -f "\$vcf_file" >> ${file_list}
    done

    if [ ! -s ${file_list} ]; then
        echo "Error: No filtered VCF files found or listed in '${file_list}' for ${meta.cohort} chromosome ${meta.chr}." >&2
        exit 1
    fi

    echo "File list for merging (${file_list}):"
    cat ${file_list}

    ${params.bcftools} concat \
        --allow-overlaps \
        --file-list ${file_list} \
        -Oz -o ${out_vcf}

    echo "Indexing merged file: ${out_vcf}"
    ${params.tabix} -p vcf ${out_vcf}
    """

    // stub:
    // // Use task.ext here too
    // """
    // #!/bin/bash
    // echo "Stub MERGE_VCF for ${task.ext.cohort} chr ${task.ext.chr}"
    // touch ${task.ext.cohort}_chr${task.ext.chr}.merged.vcf.gz
    // touch ${task.ext.cohort}_chr${task.ext.chr}.merged.vcf.gz.tbi
    // """
}
