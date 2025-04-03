// modules/phase_shapeit5.nf
// Phase genotypes using SHAPEIT5 (Phase stage)

process PHASE_SHAPEIT5 {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}/phased_shapeit5", mode: 'copy'

    // SHAPEIT5 needs bcftools for input conversion/prep if needed
    // Assuming SHAPEIT5 path is handled by beforeScript in config

    input:
    tuple val(meta), path(bcf_file), path(bcf_index) // Input is BCF from VCF_TO_BCF

    output:
    tuple val(meta), path("${meta.id}_phased_shapeit5.bcf"), path("${meta.id}_phased_shapeit5.bcf.csi"), emit: phased_bcf // SHAPEIT5 outputs BCF

    script:
    // Define output_prefix for use within the script commands
    def output_prefix = "${meta.id}_phased_shapeit5"
    // Define output_bcf based on output_prefix for clarity in script
    def output_bcf = "${output_prefix}.bcf"

    // SHAPEIT5 requires a genetic map, similar to Eagle
    def ref_map_chr_path = params.gmap.replace("CHR", meta.chr)
    def ref_map_chr = file(ref_map_chr_path)
    if (!ref_map_chr.exists()) {
        error "Genetic map file for chromosome ${meta.chr} not found at expected path: ${ref_map_chr_path}"
    }

    // SHAPEIT5 doesn't need a reference panel for phasing target samples
    // It might need chromosome length info if FASTA index isn't available? Check docs.
    // --thread is used for number of threads
    """
    echo "Starting phasing for ${meta.id} with SHAPEIT5 using ${task.cpus} threads"
    echo "SHAPEIT5 Version Info:"
    shapeit5 --version || true // Or relevant version command

    ${params.bcftools} --version

    # SHAPEIT5 Phase command
    # Requires --map, --input (target VCF/BCF), --region
    shapeit5 --phase \
        --input ${bcf_file} \
        --map ${ref_map_chr} \
        --output ${output_bcf} \
        --region ${meta.chr} \
        --thread ${task.cpus} \
        --log ${output_prefix}.log
        # Add other relevant SHAPEIT5 options

    # Validate and index output BCF
    if [ -s "${output_bcf}" ]; then
        echo "Indexing SHAPEIT5 output BCF: ${output_bcf}"
        ${params.bcftools} index -f --threads ${task.cpus} ${output_bcf}
    else
        echo "ERROR: SHAPEIT5 output file '${output_bcf}' not found or is empty!" >&2
        exit 1
    fi

    echo "SHAPEIT5 phasing completed for ${meta.id}"
    """
} 