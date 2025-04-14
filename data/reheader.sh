#!/bin/bash

for i in {1..3}; do
    cd "Cohort_test_${i}" || exit
    for j in {10..11}; do
        name="test_cohort${i}_filtered_chr${j}"
        old="${name}.recode.vcf.gz"
        new="${name}.header.recode.vcf.gz"

        # Extract the header
        bcftools view -h "$old" > "header_chr${j}.txt"

        # Modify the header.txt
        sed -i 's/^##fileformat=.*/##fileformat=VCFv4.2/' "header_chr${j}.txt"
        sed -i 's/^##FORMAT=<ID=GL,Number=.,/##FORMAT=<ID=GL,Number=G,/' "header_chr${j}.txt"

        # Apply the modified header
        bcftools reheader -h "header_chr${j}.txt" -o "$new" "$old"

        # Replace the old file with the new one
        mv "$new" "$old"

        # Index the VCF file
        tabix -p vcf "$old"
    done
    cd .. || exit
done
