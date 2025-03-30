// Process 4: CONVERT_REFP
// Convert phased VCF to M3VCF format for use as reference panel
process CONVERT_REFP {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/reference_panel", mode: 'copy'
    
    conda "bioconda::minimac4=4.1.4"

    input:
    tuple val(meta), path(phased_vcf), path(phased_vcf_index)

    output:
    tuple val(meta), path("${meta.id}_refpanel.msav"), emit: reference

    script:
    def prefix = "${meta.id}_refpanel"
    """
    echo "Converting ${phased_vcf} to Reference Panel format"
    
    # Convert to M3VCF format using minimac4
    minimac4 --compress-reference ${phased_vcf} -o ${prefix}.msav
    
    echo "Reference Panel conversion completed for ${meta.id}"
    """
} 