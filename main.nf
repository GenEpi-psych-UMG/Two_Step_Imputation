#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.help=false


// Subworkflows
include { PRE_PHASE } from './subworkflows/pre_phase'
include { PHASE_WORKFLOW } from './subworkflows/phase_workflow'
include { IMPUTATION } from './subworkflows/imputation_1st'
include { POST_IMPUTATION } from './subworkflows/post_imputation'

def printHeader() {
    log.info """
    ╔═══════════════════════════════════════════════════════════════════╗
    ║  ████████╗██╗    ██╗ ██████╗     ███████╗████████╗███████╗██████╗ ║
    ║  ╚══██╔══╝██║    ██║██╔═══██╗    ██╔════╝╚══██╔══╝██╔════╝██╔══██╗║
    ║     ██║   ██║ █╗ ██║██║   ██║    ███████╗   ██║   █████╗  ██████╔╝║
    ║     ██║   ██║███╗██║██║   ██║    ╚════██║   ██║   ██╔══╝  ██╔═══╝ ║
    ║     ██║   ╚███╔███╔╝╚██████╔╝    ███████║   ██║   ███████╗██║     ║
    ║     ╚═╝    ╚══╝╚══╝  ╚═════╝     ╚══════╝   ╚═╝   ╚══════╝╚═╝     ║
    ║                                                                   ║
    ║             IMPUTATION PIPELINE v${workflow.manifest.version}                              ║
    ║                   POWERED BY NEXTFLOW                             ║
    ╠═══════════════════════════════════════════════════════════════════╣
    ║  Input/Output Parameters:                                         ║
    ║  • Cohorts CSV   : ${params.cohorts_csv}
    ║  • Output Dir    : ${params.outdir}
    ║                                                                   ║
    ║  Reference Data:                                                  ║
    ║  • Reference FASTA: ${params.fasta}
    ║                                                                   ║
    ║  Execution Parameters:                                            ║
    ║  • Max CPUs           : ${params.max_cpus}
    ║  • Task CPUs          : ${params.max_cpus}
    ║  • Max memory         : ${params.max_cpus}
    ║  • Task Memory        : ${params.max_cpus}
    ║  • Phasing engine     : ${params.max_cpus}
    ║  • Imputation engine  : ${params.imputer}
    ╠═══════════════════════════════════════════════════════════════════╣
    ║  PIPELINE DESCRIPTION:                                            ║
    ║  Cross-imputation between cohorts where each cohort is imputed    ║
    ║  using all other cohorts as reference panels. Reference panels    ║
    ║  are created dynamically during pipeline execution.               ║
    ╚═══════════════════════════════════════════════════════════════════╝
    Pipeline started: ${new Date().format('dd-MM-yyyy HH:mm:ss')}
    """.stripIndent()
}

def printHelp() {
    log.info """
    Usage:

    nextflow run <your_repo/main.nf> -profile <standard/test/cluster> [options]

    Required Parameters:
    --cohorts_csv   Path to CSV file describing cohorts (see README).
    --outdir        Directory where results will be saved.
    --fasta         Path to reference genome FASTA file.
    --gmap_eagle2   Path to genetic map file (for Eagle/SHAPEIT5/Beagle).

    Optional Parameters:
    --phaser        Phasing tool to use ('eagle' or 'shapeit5') (Default: ${params.phaser}).
    --imputer       Imputation tool to use ('minimac4' or 'beagle') (Default: ${params.imputer}).
    --snp_select_r  Path to SNP_Selection.R script (Default: ${params.snp_select_r}).
    --extract_info_pl Path to extract_info_VCF.pl script (Default: ${params.extract_info_pl}).
    --cpus          Base number of CPUs for tasks (Default: ${params.cpus}).
    --minimac4_rounds        Number of phasing rounds for Minimac4 (Default: ${params.minimac4_rounds}).
    --chromosomes   Range of chromosomes to process (e.g., \"1..22\") (Default: 1..22).

    Tool Path Parameters (Required if not in PATH/handled by Conda):
    --eagle2        Path to Eagle executable.
    --shapeit5      Path to SHAPEIT5 executable (phase_common).
    --minimac4      Path to Minimac4 executable.
    --beagle5       Path to Beagle JAR file.
    --bref3         Path to BREF3 JAR file.
    --bcftools      Path to bcftools executable.
    --tabix         Path to tabix executable.
    --bgzip         Path to bgzip executable.

    Other Options:
    -profile        Configuration profile to use (e.g., standard, test).
    -resume         Resume previous run from cache.
    --help          Show this help message.
    """ .stripIndent()
}


workflow  {
    printHeader()

    if (params.help) {
        printHelp()
        exit 0
    }

    // Create a channel for ordered cohorts first
    ch_ordered_cohorts = Channel.fromPath(params.cohorts_csv)
                                .splitCsv(header:true, sep:',')
                                .map { row -> row.Cohort }
                                .toList()
                                .map { cohorts -> cohorts.unique() }

    // Create channel for input VCFs
    ch_inputVCFs = Channel.fromPath(params.cohorts_csv)
                                .splitCsv(header:true, sep:',')
                                .map { row ->
                                    def cohort = row.Cohort
                                    def path = row.Path.trim()
                                    def prefix = row.Prefix.trim()
                                    def suffix = row.Suffix.trim()
                                    
                                    def chromosomes = params.chromosomes ?: (1..22)
                                    
                                    return chromosomes.collect { chr ->
                                            def meta = [
                                                id: "${cohort}_chr${chr}",
                                                cohort: cohort,
                                                chr: chr
                                            ]
                                            def vcf_path = "${path}/${prefix}${chr}${suffix}"
                                            def vcf_file = file(vcf_path)
                                            if (vcf_file.exists()) {
                                                log.info "Found input VCF: ${vcf_file}"
                                                return [meta, vcf_file]
                                            } else {
                                                return null
                                            } 
                                            }
                                }
                                .flatMap()
                                .filter { it != null }
                                .ifEmpty { exit 1, "No input VCF found based on ${params.cohorts_csv}. Give me a decent CSV!! (e.g., path/prefixCHRsuffix)." }
    
    ch_bcf = PRE_PHASE(ch_inputVCFs)
    ch_phasedVCFs = PHASE_WORKFLOW(ch_bcf)
    ch_imputedVCFs = IMPUTATION(ch_phasedVCFs)
    ch_mergedVCFs = POST_IMPUTATION(ch_imputedVCFs, ch_ordered_cohorts)
    
    workflow.onError {
        log.error "Pipeline Failed!"
        if (workflow.errorMessage) {
            log.error "Error message: ${workflow.errorMessage}"
        }
        if (workflow.exitStatus) {
            log.error "Exit status: ${workflow.exitStatus}"
        }
    }
    
}