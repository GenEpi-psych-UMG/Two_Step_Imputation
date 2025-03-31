#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include modules
include { ALIGNMENT } from './modules/alignment'
include { VCF_TO_BCF } from './modules/vcf_to_bcf'
include { PHASE } from './modules/phase'
include { CONVERT_REFP } from './modules/convert_refp'
include { IMPUTE } from './modules/impute'
include { PREPARE_R_INPUT } from './modules/prepare_r_input'
include { RUN_R_SELECT } from './modules/run_r_select'
include { FILTER_VCF } from './modules/filter_vcf'
include { MERGE_VCF } from './modules/merge_vcf'

// --- Log Pipeline Header ---
def printHeader() {
    log.info """
         Two Step Imputation Pipeline - Nextflow
         ===========================================
         Cohorts CSV   : ${params.cohorts_csv}
         Output Dir    : ${params.outdir}
         Ref FASTA     : ${params.fasta}
         Genetic Map   : ${params.gmap}
         CPUs          : ${params.cpus}
         Minimac Rounds: ${params.rounds}
         ===========================================
         NOTE: This pipeline performs cross-imputation 
         between cohorts. Each cohort is imputed using 
         all other cohorts as reference panels.
         Reference panels are created during the pipeline run.
         ===========================================
         """ .stripIndent()
}

// --- Workflow Definition ---
workflow {

    printHeader()
    
    // --- Validate Parameters ---
    // Ensure required reference files exist
    if (!params.fasta || !file(params.fasta).exists()) {
        exit 1, "Reference FASTA file not found: ${params.fasta ?: 'parameter not set'}"
    }
    if (!params.gmap || !file(params.gmap).exists()) {
        exit 1, "Genetic Map file not found: ${params.gmap ?: 'parameter not set'}"
    }

    // --- Input Channel Creation ---
    // Channel to read Cohorts_Info.csv and find input VCF files per chromosome
    ch_input_vcfs = Channel.fromPath(params.cohorts_csv)
        .splitCsv(header:true, sep:',') // each row is a list of named columns
        .map { row ->
            def cohort = row.Cohort
            def pathway = row.Pathway.trim()
            def prefix = row.Prefix.trim()
            def suffix = row.Suffix.trim()
            // Generate chromosome list (1-22)
            def chromosomes = (10..11)
            // Create [meta, vcf_file] tuples for each existing VCF
            // collect operator creates a list of tuples for each chromosome in the initiated list
            return chromosomes.collect { chr ->
                def meta = [
                    id: "${cohort}_chr${chr}",
                    cohort: cohort,
                    chr: chr
                ]
                def vcf_path = "${pathway}/${prefix}${chr}${suffix}"
                def vcf_file = file(vcf_path) // Create file object
                if (vcf_file.exists()) {
                    log.info "Found input VCF: ${vcf_file}"
                    return [ meta, vcf_file ]
                } else {
                    return null
                }
            }
        }
        // output nested list of tuples created by the collect operator
        //    [
        //    [ tuple1, tuple2, ..., tuple22 ],  // For cohort1
        //    [ tuple1, tuple2, ..., tuple22 ],  // For cohort2
        //    [ tuple1, tuple2, ..., tuple22 ]   // For cohort3
        //    ]
        .flatMap() // Flatten the list of lists into individual [meta, vcf_path] emissions
        .filter { it != null } // Remove null entries where files didn't exist
        .ifEmpty { exit 1, "No input VCF files found based on ${params.cohorts_csv}. Please check paths and naming convention (e.g., path/prefixCHRsuffix)." }

    // Pipeline steps
    ALIGNMENT(ch_input_vcfs)
    VCF_TO_BCF(ALIGNMENT.out.vcf)
    PHASE(VCF_TO_BCF.out.bcf)
    CONVERT_REFP(PHASE.out.phased_vcf)
    
    // Create imputation cohort combinations
    ch_cohorts = ch_input_vcfs
        .map { meta, _vcf -> meta.cohort }
        .unique()
        .collect()
        .map { cohorts -> 
            def pairs = []
            cohorts.eachWithIndex { source, i ->
                cohorts.eachWithIndex { target, j ->
                    if (i != j) {
                        pairs.add([source, target])
                    }
                }
            }
            return pairs
        }
        .flatMap()
        
    // Combine phased VCFs with target cohorts for imputation
    ch_impute_input = PHASE.out.phased_vcf
        .combine(ch_cohorts)
        .filter { meta, _vcf, _vcf_index, source, _target ->
            return meta.cohort == source
        }
        .map { meta, vcf, vcf_index, _source, target ->
            def meta_clone = meta.clone()
            meta_clone.target_cohort = target
            return [meta_clone, vcf, vcf_index]
        }
    
    // Get reference panel files per chromosome and target cohort
    ch_ref_panels = CONVERT_REFP.out.reference
        .map { meta, msav -> 
            // Create a string key like "10_Cohort_test_1" for easy joining
            def key = "${meta.chr}_${meta.cohort}"
            return [key, msav]
        }

    // Cross-reference the reference panels with imputation inputs using join
    ch_impute_with_ref = ch_impute_input
        .map { meta, vcf, vcf_index -> 
            // Create the same format key based on target cohort
            def key = "${meta.chr}_${meta.target_cohort}"
            return [key, meta, vcf, vcf_index]
        }
        .join(ch_ref_panels, failOnMismatch: false) // Join by key, allowing mismatches
        .filter { it.size() == 5 } // Filter out any mismatches (join will return fewer elements)
        .map { key, meta, vcf, vcf_index, ref_panel ->
            return [meta, vcf, vcf_index, ref_panel]
        }
    
    // Run imputation with matched reference panels
    IMPUTE(ch_impute_with_ref)
    
    // Extract info from imputed VCF files
    PREPARE_R_INPUT(IMPUTE.out.imputed_vcf)
    
    // Group info files by cohort for R analysis
    ch_r_input = PREPARE_R_INPUT.out.info_files
        .map { meta, info_file -> 
            return [meta.cohort, meta, info_file]
        }
        .groupTuple(by: 0)
        .map { _cohort, _metaItems, infoFiles -> 
            // Create per-cohort list file with paths to all info files
            def meta = [
                id: "${_cohort}_R_input",
                cohort: _cohort
            ]
            return [meta, infoFiles]
        }
    
    // Run R script to select SNPs
    RUN_R_SELECT(ch_r_input)
    
    // Prepare for VCF filtering
    ch_filter_input = IMPUTE.out.imputed_vcf
        .filter { meta, _vcf -> meta != null && meta.cohort != null }           // <-- New filter ensuring meta is defined
        .map { meta, vcf -> [ meta.cohort, meta, vcf ] }
        .combine( RUN_R_SELECT.out.filter_lists
            .filter { meta, _list -> meta != null && meta.cohort != null }     // <-- New filter in R channel
            .map { meta, list -> [ meta.cohort, list ] }, by: 0)
        .map { _cohort, meta, vcf, filter_list -> [ meta, vcf, filter_list ] }
    
    // Filter VCF files
    FILTER_VCF(ch_filter_input)
    
    // Group filtered VCFs by cohort and chromosome for merging
    ch_merge_input = FILTER_VCF.out.filtered_vcf
        .filter { meta, _vcf -> meta != null && meta.cohort != null && meta.chr != null }  // <-- New filter to skip null meta
        .map { meta, vcf -> 
            return [ meta.cohort, meta.chr, meta, vcf ]
        }
        .groupTuple(by: [0, 1]) // Group by [cohort, chr]
        .map { _cohort, _chr, _metas, vcfs ->
            def meta = [
                id: "${_cohort}_chr${_chr}_merged",
                cohort: _cohort,
                chr: _chr
            ]
            return [ meta, vcfs ]
        }
    
    // Merge filtered VCFs from different reference cohorts for each chromosome
    MERGE_VCF(ch_merge_input)
}