#!/bin/bash

# Source the configuration file to load paths
source "$1"

# Check if the config.sh file was sourced properly
if [ -z "$eagle" ] || [ -z "$gmap" ] || [ -z "$cpus" ] || [ -z "$rounds" ]; then
  echo "Error: Config file variables not loaded properly."
  exit 1
fi

# Check if the config.sh file was sourced properly
if [ -z "$eagle" ] || [ -z "$gmap" ] || [ -z "$cpus" ] || [ -z "$rounds" ]; then
  echo "Error: parameters needed for this script are not found.."
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
startalltime=`date +%s`
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


for chr in {1..23}; do
starttime=`date +%s`
input_bcf="${pathway}/${Cohort_prefix}_chr${chr}_aligned.bcf.gz"
output_phased="${pathway}/${Cohort_prefix}_chr${chr}_phased"
output_convert="${pathway}/${Cohort_prefix}_chr${chr}_phased_RefP.msav"
 # Debugging: Print the file paths being used
    echo "Input VCF: $input_bcf"
    echo "Output VCF: $output_phased"

 # Check if the input VCF file exists before processing
    if [ ! -f "$input_bcf" ]; then
      echo "Warning: Input file '$input_bcf' not found. Skipping."
      continue
    fi


"$eagle" \
 --vcf "$input_bcf" \
 --Kpbwt=10000 \
 --pbwtIters="$rounds" \
 --geneticMapFile="$gmap" \
 --vcfOutFormat=z \
 --outPrefix "$output_phased" \
 --numThreads "$cpus" 

tabix "$output_phased".vcf.gz

# Check if bcftools was successful
    if [ $? -eq 0 ]; then
      echo "phasing for ${input_bcf} is done"
    else
      echo "Error: Phasing for ${input_bcf} failed"
    fi
endtime=`date +%s`
echo "runtime of phase.sh was $((endtime-starttime)) seconds for chr ${chr}"
  done
endalltime=`date +%s`
echo "runtime of phase.sh was $((endalltime-startalltime)) seconds for ${Cohort_prefix} "
     echo "phasing for ${Cohort_prefix} done"
done
