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
        .splitCsv(header:true, sep:',')
        .map { row ->
            def cohort = row.Cohort
            def pathway = row.Pathway.trim() // Trim potential whitespace
            def prefix = row.Prefix.trim()
            def suffix = row.Suffix.trim()
            // Generate chromosome list (1-22 and X)
            def chromosomes = (1..22) + ['X']
            // Create [meta, vcf_file] tuples for each existing VCF
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
                    log.warn "Input VCF not found, skipping: ${vcf_path}"
                    return null // Return null if file doesn't exist
                }
            }
        }
        .flatMap() // Flatten the list of lists into individual [meta, vcf_path] emissions
        .filter { it != null } // Remove null entries where files didn't exist
        .ifEmpty { exit 1, "No input VCF files found based on ${params.cohorts_csv}. Please check paths and naming convention (e.g., path/prefixCHRsuffix)." }

    // Pipeline steps
    ALIGNMENT(ch_input_vcfs)
    VCF_TO_BCF(ALIGNMENT.out.vcf)
    PHASE(VCF_TO_BCF.out.bcf)
    CONVERT_REFP(PHASE.out.phased_vcf)
    
    // Create imputation cohort combinations
    // Extract unique cohorts from the input channel
    ch_cohorts = ch_input_vcfs
        .map { meta, _vcf -> meta.cohort }
        .unique()
        .collect()
        .map { cohorts -> 
            // Generate all cohort pairs for cross-imputation (excluding self pairs)
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
        .filter { meta, _vcf, _vcf_index, cohort_pair -> 
            // Only keep pairs where the VCF cohort matches the source cohort in the pair
            return meta.cohort == cohort_pair[0]
        }
        .map { meta, vcf, vcf_index, cohort_pair -> 
            // Convert to [meta with target, vcf]
            def new_meta = meta.clone()
            new_meta.target_cohort = cohort_pair[1]
            return [new_meta, vcf, vcf_index]
        }
    
    // Get reference panel files per chromosome and target cohort
    // First, create a reusable reference panel channel
    ch_ref_panels_all = CONVERT_REFP.out.reference
        .map { meta, msav -> 
            // Key: [chromosome, cohort]
            return [[chr: meta.chr, cohort: meta.cohort], msav]
        }
        .groupTuple(by: 0)
        .map { key, files -> 
            if (files.size() > 1) {
                log.warn "Multiple reference files found for chr ${key.chr} cohort ${key.cohort}, using first: ${files[0]}"
            }
            return [key, files[0]]
        }
        .collect() // Collect all reference panels

    // Cross-reference the reference panels with imputation inputs
    ch_impute_with_ref = ch_impute_input
        .combine(ch_ref_panels_all)
        .map { meta, vcf, vcf_index, ref_panels -> 
            // Find the matching reference panel for this imputation
            def target_key = [chr: meta.chr, cohort: meta.target_cohort]
            def matched_ref_panel = null
            
            // Find the matching reference panel in the collected list
            ref_panels.each { ref_key, ref_file ->
                if (ref_key.chr == target_key.chr && ref_key.cohort == target_key.cohort) {
                    matched_ref_panel = ref_file
                }
            }
            
            if (matched_ref_panel == null) {
                log.error "No reference panel found for ${meta.id} with target cohort ${meta.target_cohort}"
                return null
            }
            
            // Return input with matching reference panel
            return [meta, vcf, vcf_index, matched_ref_panel]
        }
        .filter { it != null }
    
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
    // We need to combine each imputed VCF with its corresponding filter list
    ch_filter_input = IMPUTE.out.imputed_vcf
        .map { meta, vcf -> [meta.cohort, meta, vcf] }
        .combine(RUN_R_SELECT.out.filter_lists.map { meta, list -> [meta.cohort, list] }, by: 0)
        .map { _cohort, meta, vcf, filter_list -> [meta, vcf, filter_list] }
    
    // Filter VCF files
    FILTER_VCF(ch_filter_input)
    
    // Group filtered VCFs by cohort and chromosome for merging
    ch_merge_input = FILTER_VCF.out.filtered_vcf
        .map { meta, vcf -> 
            return [meta.cohort, meta.chr, meta, vcf]
        }
        .groupTuple(by: [0, 1]) // Group by [cohort, chr]
        .map { _cohort, _chr, _metas, vcfs ->
            def meta = [
                id: "${_cohort}_chr${_chr}_merged",
                cohort: _cohort,
                chr: _chr
            ]
            return [meta, vcfs]
        }
    
    // Merge filtered VCFs from different reference cohorts for each chromosome
    MERGE_VCF(ch_merge_input)
} 