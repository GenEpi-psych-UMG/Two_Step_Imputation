#!/bin/bash

# Source the configuration file to load paths
source "$1"

# Check if the config.sh file was sourced properly
if [ -z "$fasta" ] || [ -z "$cohorts" ]; then
  echo "Error: Config file variables not loaded properly."
  exit 1
fi

# Check if the fasta file exists
if [ ! -f "$fasta" ]; then
  echo "Error: Fasta file not found."
  exit 1
fi

# Initialize arrays for pathway, filename, cohort, and suffix
Pathways=()
Prefixes=()
Cohorts=()
Suffixes=()

# Read prefixes and pathways from the cohort info file (assuming it's CSV)
skip_header=true
while IFS=, read -r Cohort Pathway Prefix Suffix; do
if $skip_header; then
    skip_header=false
    continue
  fi
  Pathways+=("$Pathway")
  Prefixes+=("$Prefix")
  Cohorts+=("$Cohort")
  Suffixes+=("$Suffix")
done < "$cohorts"

# Iterate over prefixes and pathways for chromosomes 1 to 22
for i in "${!Cohorts[@]}"; do
  pathway=${Pathways[$i]}
  Input_prefix=${Prefixes[$i]}
  Cohort_prefix=${Cohorts[$i]}
  Suffix=${Suffixes[$i]}

# Trim any leading/trailing whitespace from variables
  pathway=$(echo "$pathway" | sed 's/[[:space:]]*$//')
  #Input_prefix=$(echo "$Input_prefix" | xargs)
  Cohort_prefix=$(echo "$Cohort_prefix" | xargs)
  #Suffix=$(echo "$Suffix" | xargs)

# Debugging: Print the values of the variables being used
  echo "Processing Cohort: $Cohort_prefix"
  echo "Pathway: $pathway"
  #echo "Input Prefix: $Input_prefix"
  #echo "Suffix: $Suffix"

  for chr in {1..22}; do
    input_vcf="${pathway}/${Cohort_prefix}_chr${chr}_aligned.vcf.gz"
    output_bcf="${pathway}/${Cohort_prefix}_chr${chr}_aligned.bcf.gz"

 # Debugging: Print the file paths being used
    echo "Input VCF: $input_vcf"
    echo "Output VCF: $output_bcf"

    # Check if the input VCF file exists before processing
    if [ ! -f "$input_vcf" ]; then
      echo "Warning: Input file '$input_vcf' not found. Skipping."
      continue
    fi

    "$Bcftools" view -c 1 -O b -o  "$output_bcf" "$input_vcf"


"$Tabix" -p vcf "$output_bcf"
echo "chromosome $chr done"

    
    # Check if bcftools was successful
    if [ $? -eq 0 ]; then
      echo "BCF Conversion for ${input_vcf} is done"
    else
      echo "Error: Conversion for ${input_vcf} failed"
    fi
  done
     echo "BCF Conversion for ${Cohort_prefix} done"
done
