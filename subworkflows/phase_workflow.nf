#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Include modules
include { EAGLE2 } from '../modules/phase_eagle2'
include { SHAPEIT5 } from '../modules/phase_shapeit5'
include { BCF_TO_VCF } from '../modules/bcf_to_vcf'
include { gmapFinder } from '../lib/gmap_finder'

workflow PHASE_WORKFLOW {
    take:
    ch_bcf // Channel of [meta, bcf_file, bcf_index]

    main:  
    ch_phasedVCFs = Channel.empty()

    // --- Phasing Step (Conditional) ---
    if (params.phaser == 'eagle2') {
        EAGLE2(ch_bcf)
        ch_phasedVCFs = EAGLE2.out.phased_vcf // Eagle outputs VCF.gz
    
    
    
    } else if (params.phaser == 'shapeit5') {
        
        ref_map_chr = gmapFinder(params.phaser, params.gmap_shapeit5)
        
        // ref_map_chr.view { "Genetic map file: ${it}" }

        ch_bcf_keyed = ch_bcf
                        .map { meta, bcf_file, bcf_index -> 
                        [meta.chr.toString(), meta, bcf_file, bcf_index] 
                        }
    
        // Join with the ref_map_chr channel using chromosome as key
        ch_bcf_keyed.combine(ref_map_chr, by: 0)
                    .map { chr, meta, bcf_file, bcf_index, ref_map_file -> 
            [meta, bcf_file, bcf_index, ref_map_file] 
            }
            // .view { "SHAPEIT5 input: ${it}" }
            .set { ch_bcf_with_map }

        
        SHAPEIT5(ch_bcf_with_map)
        // SHAPEIT5.out.view { "SHAPEIT5 output: ${it}" }
        BCF_TO_VCF(SHAPEIT5.out.phased_bcf) // Convert SHAPEIT5 BCF to VCF
        // BCF_TO_VCF.out.view { "BCF_TO_VCF output: ${it}" }
        ch_phasedVCFs = BCF_TO_VCF.out.converted_vcf
    } else {
        error "Invalid phaser specified: ${params.phaser}. Choose 'eagle' or 'shapeit5'."
    }

    // ch_phasedVCFs.view { "Phased VCF ready for Imputation/RefPanel: ${it[0].id}" }

    emit:
    phased_vcf = ch_phasedVCFs // Output [meta, vcf, index]
}