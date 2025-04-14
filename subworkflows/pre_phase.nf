#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Include modules
include { ALIGNMENT } from '../modules/alignment'
include { VCF_TO_BCF } from '../modules/vcf_to_bcf'

workflow PRE_PHASE {
    take:
    ch_inputVCFs // Channel of [meta, vcf_file]

    main:
    ALIGNMENT(ch_inputVCFs)
    VCF_TO_BCF(ALIGNMENT.out.vcf)

    emit:
    bcf = VCF_TO_BCF.out.bcf // Output [meta, bcf_file, bcf_index]
}