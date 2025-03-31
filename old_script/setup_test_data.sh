#!/bin/bash

# Setup script for Two-Step Imputation Pipeline test data
echo "Setting up test data directories for Two-Step Imputation Pipeline..."

# Create test directories if they don't exist
mkdir -p Cohort_test/Cohort_test_1
mkdir -p Cohort_test/Cohort_test_2
mkdir -p Cohort_test/Cohort_test_3

# Define the chromosomes to create dummy data for
chrs=("1" "2" "3" "X")

# Create dummy VCF files for each test cohort and chromosome
for cohort in {1..3}; do
  for chr in "${chrs[@]}"; do
    dummy_vcf="Cohort_test/Cohort_test_${cohort}/chr${chr}.vcf.gz"
    
    # Only create if it doesn't exist
    if [ ! -f "$dummy_vcf" ]; then
      echo "Creating dummy VCF: $dummy_vcf"
      
      # Create a minimal dummy VCF file
      echo -e '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1' > temp.vcf
      
      # Add 10 dummy variants
      for pos in {1000..10000..1000}; do
        echo -e "chr${chr}\t${pos}\trs${pos}\tA\tT\t100\tPASS\t.\tGT\t0/1" >> temp.vcf
      done
      
      # Compress and index
      bgzip -c temp.vcf > "$dummy_vcf"
      tabix -p vcf "$dummy_vcf"
      rm temp.vcf
    fi
  done
done

echo "Test data setup complete!"
echo ""
echo "To run the pipeline with test data:"
echo "1. Update the reference genome path in nextflow.config:"
echo "   fasta = '/path/to/human_g1k_v37.fasta'"
echo "   gmap = '/path/to/genetic_map_hg19.txt.gz'"
echo ""
echo "2. Run the pipeline:"
echo "   nextflow run main.nf"
echo ""
echo "Results will be stored in: ../Two_Step_Imputation_Results" 