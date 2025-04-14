#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Include modules
include { CONVERT_REFP as CONVERT_REFP_MSAV } from '../modules/convert_refp'
include { CONVERT_REFP as CONVERT_REFP_BREF3 } from '../modules/convert_refp_bref3'
include { pairImputeJobs} from '../lib/impute_job_builder'
include { gmapFinder } from '../lib/gmap_finder'
include { MINIMAC4 } from '../modules/impute_minimac4'
include { BEAGLE5 } from '../modules/impute_beagle5'

workflow IMPUTATION {
    take:
    ch_phasedVCFs // [meta, vcf, index] from phasing
    
    main:

    
    if (params.imputer == "minimac4") {
        CONVERT_REFP_MSAV(ch_phasedVCFs)
        
        pairImputeJobs(ch_phasedVCFs, CONVERT_REFP_MSAV.out.reference).set { ch_impute_jobs }
        
        MINIMAC4(ch_impute_jobs)
        ch_imputedVCFs = MINIMAC4.out.imputed_vcf
    } 
    else if (params.imputer == "beagle5") {
        ch_java_memory = (params.memory =~ /(\d+)(?:\.?\s*GB)/)[0][1].toInteger()
        
        // For Beagle imputation
        CONVERT_REFP_BREF3(ch_phasedVCFs, ch_java_memory)

        pairImputeJobs(ch_phasedVCFs, CONVERT_REFP_BREF3.out.reference)
            .map { newMeta, phased_vcf, phased_index, ref_panel -> 
            [newMeta.chr.toString(), newMeta, phased_vcf, phased_index, ref_panel]
            }
            .set { ch_impute_jobs_keyed }

        ref_map_chr = gmapFinder(params.imputer, params.gmap_beagle5)

        ch_impute_jobs_keyed.combine(ref_map_chr, by: 0)
            .map { chr, meta, phased_vcf, phased_index, ref_panel, ref_map_file -> 
                [meta, phased_vcf, phased_index, ref_panel, ref_map_file]
            }
            .view { "BEAGLE5 input: ${it}" }
            .set { ch_impute_jobs }
        
        BEAGLE5(ch_impute_jobs, ch_java_memory)
        
        // Ensure the output has the same structure as IMPUTE_MINIMAC4
        ch_imputedVCFs = BEAGLE5.out.imputed_vcf

    }
    else {
        error "Unsupported imputer: ${params.imputer}. Choose 'minimac4' or 'beagle5'."
    }

    emit:
    imputed_vcf = ch_imputedVCFs // [meta, vcf, index] format for POST_IMPUTATION
}