// Process 2: VCF_TO_BCF
// Convert VCF to BCF format for faster processing
process VCF_TO_BCF {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/bcf", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("${meta.id}.bcf.gz"), path("${meta.id}.bcf.gz.tbi"), emit: bcf

    script:
    def prefix = "${meta.id}"
    """
    echo "Converting ${vcf} to BCF format"
    
    # Convert to BCF with call filter (-c 1)
    bcftools view -c 1 -O b -o ${prefix}.bcf.gz ${vcf}
    
    # Index the output
    bcftools -t ${prefix}.bcf.gz
    
    echo "BCF conversion completed for ${meta.id}"
    """
} 