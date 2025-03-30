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
         Eagle Path    : ${params.eagle}
         Ref Panel VCF : ${params.ref_panel_m3vcf ?: 'Not Provided - REQUIRED for Imputation'}
         CPUs          : ${params.cpus}
         Minimac Rounds: ${params.rounds}
         ===========================================
         """ .stripIndent()
}

// --- Workflow Definition ---
workflow {

    printHeader()
    
    // --- Validate Parameters ---
    // Ensure required reference files and reference panel are specified
    if (!params.fasta || !file(params.fasta).exists()) {
        exit 1, "Reference FASTA file not found: ${params.fasta ?: 'parameter not set'}"
    }
    if (!params.gmap || !file(params.gmap).exists()) {
        exit 1, "Genetic Map file not found: ${params.gmap ?: 'parameter not set'}"
    }
    if (!params.ref_panel_m3vcf) {
        exit 1, "Parameter 'params.ref_panel_m3vcf' (path/glob to reference panel m3vcf files) is required."
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

    // --- Reference Panel Channel ---
    // Create a channel mapping chromosome to its reference panel file
    ch_ref_panel_map = Channel.fromPath(params.ref_panel_m3vcf)
        .map { file ->
            // Extract chromosome from filename (handles chr1, chr22, chrX, chrY)
            def matcher = file.getName() =~ /chr([0-9]+|[XY])/ 
            if (!matcher) {
                log.warn "Could not determine chromosome for reference panel file: ${file}. Expecting format like '...chr<NUM|X|Y>...'. Skipping."
                return null
            }
            def chr = matcher[0][1]
            log.info "Found reference panel VCF: ${file} for chr: ${chr}" 
            return [ chr, file ] // Emit tuple [ chr, path ]
        }
        .filter { it != null } 
        .groupTuple() // Group by chromosome
        .map { chr, files -> 
            if (files.size() > 1) {
                log.warn "Multiple reference files found for chromosome ${chr}, using first one: ${files[0]}"
            }
            return [ chr, files[0] ] 
        }
        .toMap() // Collect into a map [chr: path]
        .ifEmpty { exit 1, "No reference panel m3vcf files found matching glob: ${params.ref_panel_m3vcf}" }

    // For now, just view the inputs to verify the channels are working
    ch_input_vcfs.view{ meta, vcf -> "Input VCF: ${meta.id}, ${vcf}" }
    ch_ref_panel_map.view{ chr, file -> "Reference Panel: chr${chr}, ${file}" }

    // Pipeline steps
    ALIGNMENT(ch_input_vcfs)
    VCF_TO_BCF(ALIGNMENT.out.vcf)
    PHASE(VCF_TO_BCF.out.bcf)
    CONVERT_REFP(PHASE.out.phased_vcf)
    
    // Create imputation cohort combinations
    // Each cohort needs to be imputed using every other cohort as reference
    // First, get unique cohorts from the input channel
    ch_source_cohorts = ch_input_vcfs
        .map { meta, vcf -> meta.cohort }
        .unique()
        .collect()
    
    ch_target_cohorts = ch_source_cohorts.flatMap { it }
    
    // Generate cross-combinations of all cohorts without self-comparisons
    all_cohort_pairs = ch_source_cohorts
        .combine(ch_target_cohorts)
        .filter { source_cohorts, target_cohort -> 
            // Exclude self-comparisons
            return source_cohorts.contains(target_cohort) && source_cohorts.size() > 1
        }
        .flatMap { source_cohorts, target_cohort ->
            // Generate all valid pairs [source, target] excluding [target, target]
            return source_cohorts.findAll { it != target_cohort }.collect { [it, target_cohort] }
        }
    
    // Combine phased VCFs with target cohorts for imputation
    ch_impute_input = PHASE.out.phased_vcf
        .combine(all_cohort_pairs)
        .filter { meta, vcf, vcf_index, cohort_pair -> 
            // Only keep pairs where the VCF cohort matches the first cohort in the pair
            return meta.cohort == cohort_pair[0]
        }
        .map { meta, vcf, vcf_index, cohort_pair -> 
            // Convert to [meta with target, vcf]
            def new_meta = meta.clone()
            new_meta.target_cohort = cohort_pair[1]
            return [new_meta, vcf, vcf_index]
        }
    
    // Get reference panel files per chromosome and target cohort
    ch_ref_panels = CONVERT_REFP.out.reference
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
        .toMap() // Map with [chr:cohort] as key
    
    // Combine input VCFs with corresponding reference panels
    IMPUTE(ch_impute_input, ch_ref_panels)
    
    // Extract info from imputed VCF files
    PREPARE_R_INPUT(IMPUTE.out.imputed_vcf)
    
    // Group info files by cohort for R analysis
    ch_r_input = PREPARE_R_INPUT.out.info_files
        .map { meta, info_file -> 
            return [meta.cohort, [meta, info_file]]
        }
        .groupTuple()
        .map { cohort, files -> 
            // Create per-cohort list file with paths to all info files
            def metaItems = files.collect { it[0] }
            def infoFiles = files.collect { it[1] }
            def meta = [
                id: "${cohort}_R_input",
                cohort: cohort
            ]
            return [meta, infoFiles]
        }
    
    // Run R script to select SNPs
    RUN_R_SELECT(ch_r_input)
    
    // Prepare for VCF filtering
    // We need to combine each imputed VCF with its corresponding filter list
    ch_filter_input = IMPUTE.out.imputed_vcf
        .combine(RUN_R_SELECT.out.filter_lists, by: 0) // Match by cohort
    
    // Filter VCF files
    FILTER_VCF(ch_filter_input)
    
    // Group filtered VCFs by cohort and chromosome for merging
    ch_merge_input = FILTER_VCF.out.filtered_vcf
        .map { meta, vcf -> 
            return [meta.cohort, meta.chr, vcf]
        }
        .groupTuple(by: [0, 1]) // Group by [cohort, chr]
        .map { cohort, chr, vcfs ->
            def meta = [
                id: "${cohort}_chr${chr}_merged",
                cohort: cohort,
                chr: chr
            ]
            return [meta, vcfs]
        }
    
    // Merge filtered VCFs from different reference cohorts for each chromosome
    MERGE_VCF(ch_merge_input)
} 