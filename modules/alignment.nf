// Process 1: ALIGNMENT
// Aligns VCF files to reference genome
process ALIGNMENT {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_aligned.vcf.gz"), path("${meta.id}_aligned.vcf.gz.tbi"), emit: vcf

    script:
    def prefix = "${meta.id}_aligned"
    """
    echo "Running Alignment for ${meta.id}"
    
    # Normalize and align VCF against reference
    ${params.bcftools} norm -c s -f ${params.fasta} -O z -o ${prefix}.vcf.gz ${vcf}
    
    # Index the output
    ${params.tabix} -p vcf ${prefix}.vcf.gz
    
    echo "Alignment completed for ${meta.id}"
    """
} 