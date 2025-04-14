// Process 4: CONVERT_REFP
// Convert phased VCF to M3VCF format for use as reference panel
process CONVERT_REFP {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/reference_panel", mode: 'copy'

    input:
    tuple val(meta), path(phased_vcf), path(phased_vcf_index)
    val(ch_java_memory) // Memory for java task

    output:
    tuple val(meta), path("${meta.id}_refpanel.bref3"), emit: reference

    script:
    def prefix = "${meta.id}_refpanel"
    def output = "${prefix}.bref3"
    """
    echo "Converting ${phased_vcf} to Reference Panel format (bref3)"
    
    # Convert to BREF3 format using bref3 converter beagle5.5
    java -jar  -Xmx${ch_java_memory}G ${params.bref3} ${phased_vcf} > ${output}
    
    echo "Reference Panel conversion completed for ${meta.id}"
    """
} 