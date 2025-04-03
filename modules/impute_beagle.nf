// modules/impute_beagle.nf
// Impute genotypes using Beagle 5.5

process IMPUTE_BEAGLE {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/imputed_beagle/${meta.target_cohort}", mode: 'copy'


    input:
    // Beagle takes phased target VCF and reference VCF/BCF panel
    tuple val(meta), path(phased_vcf_or_bcf), path(phased_index), path(ref_panel_vcf_or_bcf) // Note: Assumes CONVERT_REFP provides VCF/BCF

    output:
    // Construct filenames using meta within the output block
    tuple val(meta), path("${meta.id}_imputed_beagle_${meta.target_cohort}.vcf.gz"), path("${meta.id}_imputed_beagle_${meta.target_cohort}.vcf.gz.tbi"), emit: imputed_vcf
    tuple val(meta), path("${meta.id}_imputed_beagle_${meta.target_cohort}.log"), emit: imputation_log // Beagle log file

    script:
    // Define output_prefix for use within the script commands
    def output_prefix = "${meta.id}_imputed_beagle_${meta.target_cohort}"

    // Beagle needs a genetic map
    def ref_map_chr_path = params.gmap.replace("CHR", meta.chr)
    def ref_map_chr = file(ref_map_chr_path)
    if (!ref_map_chr.exists()) {
        error "Genetic map file for chromosome ${meta.chr} not found at expected path: ${ref_map_chr_path}"
    }

    // Calculate appropriate memory for Java VM (e.g., 80% of task memory)
    def java_mem = task.memory ? (task.memory.toMega() * 0.8).intValue() : 4096 // Default 4GB if task memory isn't set
    if (!params.beagle) {
        error "Path to Beagle jar file not specified in params.beagle"
    }

    """
    echo "Starting imputation for ${meta.id} (Target: ${meta.target_cohort}) with Beagle using ${task.cpus} threads"
    echo "Reference panel file: ${ref_panel_vcf_or_bcf}"
    echo "Input VCF/BCF: ${phased_vcf_or_bcf}"
    echo "Java Version:"
    ${params.java} -version || true
    echo "Tabix Version:"
    tabix --version || true

    # Run Beagle imputation
    ${params.java} -Xmx${java_mem}m -jar ${params.beagle} \
        gt=${phased_vcf_or_bcf} \
        ref=${ref_panel_vcf_or_bcf} \
        map=${ref_map_chr} \
        out=${output_prefix} \
        chrom=${meta.chr} \
        impute=true \
        nthreads=${task.cpus} \
        gp=true \
        ap=true \
        seed=-99999 \
        # Add other relevant Beagle options, e.g., ne=, err=

    # Check if Beagle output VCF exists and is not empty
    if [ -s "${output_prefix}.vcf.gz" ]; then
        echo "Indexing Beagle output VCF: ${output_prefix}.vcf.gz"
        ${params.tabix} -p vcf ${output_prefix}.vcf.gz
    else
        echo "ERROR: Beagle output file '${output_prefix}.vcf.gz' not found or is empty!" >&2
        exit 1
    fi

    echo "Beagle imputation completed for ${meta.id}"
    """
} 