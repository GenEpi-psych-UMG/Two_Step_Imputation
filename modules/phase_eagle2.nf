// Process 3: PHASE
// Phase haplotypes using Eagle
process EAGLE2 {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/phased_eagle2", mode: 'copy'
    
    input:
    tuple val(meta), path(bcf), path(bcf_index)

    output:
    tuple val(meta), path("${meta.id}_phased_eagle2.vcf.gz"), path("${meta.id}_phased_eagle2.vcf.gz.tbi"), emit: phased_vcf

    script:
    def prefix = "${meta.id}_phased_eagle2"
    """
    echo "Phasing ${bcf} with Eagle"
    
    # Run Eagle for phasing
    ${params.eagle2} \
      --vcf ${bcf} \
      --Kpbwt=10000 \
      --pbwtIters=${params.eagle2_rounds} \
      --geneticMapFile=${params.gmap_eagle2} \
      --vcfOutFormat=z \
      --outPrefix ${prefix} \
      --numThreads ${task.cpus}
    
    # Index the phased output
    ${params.tabix} -p vcf ${prefix}.vcf.gz
    
    echo "Phasing completed for ${meta.id}"
    """
} 