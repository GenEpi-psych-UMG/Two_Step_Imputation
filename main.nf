#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Log Pipeline Header ---
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
    // .view { meta, vcf -> "Input Channel: ID=${meta.id}, File=${vcf}" }


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
    .toMap() // Collect into a map [chr: path]
    .ifEmpty { exit 1, "No reference panel m3vcf files found matching glob: ${params.ref_panel_m3vcf}" }

// --- Workflow Definition ---
workflow {

    main:
    // Pass input VCFs through the pipeline steps
    // We will uncomment and implement these processes one by one
    // ALIGNMENT(ch_input_vcfs)
    // VCF_TO_BCF(ALIGNMENT.out.vcf)
    // PHASE(VCF_TO_BCF.out.bcf_indexed)
    // IMPUTE(PHASE.out.phased_vcf.combine(ch_ref_panel_map, by: 0)) // Combine by chromosome key
    // PREPARE_R_INPUT(IMPUTE.out.imputed_vcf)
    // RUN_R_SELECT(PREPARE_R_INPUT.out.r_input_files)
    // FILTER_VCF(IMPUTE.out.imputed_vcf.combine(RUN_R_SELECT.out.snp_list, by: 0))
    // MERGE_VCF(FILTER_VCF.out.filtered_vcf.map { meta, vcf -> [meta.cohort, vcf] }.groupTuple())

    // For now, just view the inputs
    ch_input_vcfs.view()
    ch_ref_panel_map.view()

    // emit:
    // Define final outputs of the workflow here later
    // merged_vcfs = MERGE_VCF.out

}

// --- Process Definitions --- 
// We will add process definitions below this line

/*
process ALIGNMENT {
    tag "$meta.id Alignment"
    publishDir "${params.outdir}/${meta.cohort}/alignment", mode: 'copy', pattern: "${meta.id}_aligned.vcf.gz"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_aligned.vcf.gz"), emit: vcf

    script:
    def prefix = "${meta.id}_aligned"
    // Translate commands from Alignment.sh here
    // This likely involves bcftools norm or other alignment checks
    // Placeholder: Just copy and rename for now
    """
    echo "Running Alignment Check/Normalization for ${meta.id}"
    
    # Example using bcftools norm (Requires params.fasta to be indexed, e.g., with samtools faidx)
    # ${params.bcftools} norm -m -any -f ${params.fasta} ${vcf} -O z -o ${prefix}.vcf.gz
    
    # Placeholder: copy input if Alignment.sh logic is complex or TBD
    cp ${vcf} ${prefix}.vcf 
    gzip ${prefix}.vcf
    
    echo "Alignment step finished for ${meta.id}"
    """
}
*/

// ... Other process definitions will follow ... 