#!/bin/bash

# Source the configuration file to load paths and parameters
source "$1"

# Check if the necessary variables are loaded properly
if [ -z "$Bcftools" ] || [ -z "$cohorts" ]; then
  echo "Error: Config file variables not loaded properly."
  exit 1
fi

# Check if the cohorts file exists
if [ ! -f "$cohorts" ]; then
  echo "Error: Cohorts file not found."
  exit 1
fi

# Read the cohorts file into arrays
cohort_names=()
cohort_paths=()
Prefixes=()
Suffixes=()


skip_header=true
while IFS=',' read -r Cohort Pathway Prefix Suffix; do
if $skip_header; then
    skip_header=false
    continue
  fi
  cohort_names+=("$Cohort")
  cohort_paths+=("$Pathway")
  Prefixes+=("$Prefix")
  Suffixes+=("$Suffix")
done < "$cohorts"

# Perform pairwise imputation
for i in "${!cohort_names[@]}"; do
  input_cohort="${cohort_names[$i]}"
  input_pathway="${cohort_paths[$i]}"

  #for j in "${!cohort_names[@]}"; do
    #target_cohort="${cohort_names[$j]}"
    #target_pathway="${cohort_paths[$j]}"

# Trim any leading/trailing whitespace from variables
  input_pathway=$(echo "$input_pathway" | sed 's/[[:space:]]*$//')
  #target_pathway=$(echo "$target_pathway" | sed 's/[[:space:]]*$//')
  input_cohort=$(echo "$input_cohort" | xargs)
  #target_cohort=$(echo "$target_cohort" | xargs)

    # Ensure input cohort and target cohort are not the same
    #if [ "$input_cohort" != "$target_cohort" ]; then
mkdir -p "${input_pathway}/Outcome"
VCF_DIR="${input_pathway}"
      # Loop through chromosomes 1 to 2
      for chr in {1..22}; do

        # Define the input, target, and output file paths
        #Input_file="${input_pathway}/${input_cohort}_chr${chr}_imputed_${target_cohort}.dose.vcf.gz"
        #filter_file="${input_pathway}/Keep_list_main_chr${chr}_${target_cohort}.txt"
        Output_file="${input_pathway}/Outcome/${input_cohort}_chr${chr}_imputed_combined.vcf.gz"
        

VCF_FILES=$(find "$VCF_DIR" -name "*chr${i}*addon.recode.vcf.gz" -o -name "*chr${i}*kept.recode.vcf.gz")
        # Run the minimac4 command
#echo "Input_file= '$Input_file'"
#echo "filter_file= '$filter_file'"
echo "Output_file= '$Output_file'"

 # Change to the working directory
        if cd "${VCF_DIR}"; then
          echo "Changed to working directory: ${WD}"
        else
          echo "Error: Could not change to working directory: ${WD}"
          exit 1
        fi

VCF_FILES=$(find "$VCF_DIR" -name "*chr${chr}*addon.recode.vcf.gz" -o -name "*chr${chr}*kept.recode.vcf.gz")
if [ -z "$VCF_FILES" ]; then
    echo "No VCF files found for chromosome ${chr}. Skipping..."
    continue
  fi

  # Print the VCF files found
  echo "VCF files found for chromosome ${chr}:"
  echo "$VCF_FILES"

# Convert the VCF_FILES into a space-separated list for bcftools
  VCF_FILES_ARRAY=($VCF_FILES)
starttime=`date +%s`

  # Check if there are at least two files to concatenate
  if [ ${#VCF_FILES_ARRAY[@]} -lt 2 ]; then
    echo "Not enough VCF files to concatenate for chromosome ${chr}. Skipping..."
    continue
  fi
bcftools concat -a ${VCF_FILES_ARRAY[@]} -Oz -o "$Output_file"

  # Check if bcftools concat succeeded
  if [ $? -ne 0 ]; then
    echo "Error: Concatenation failed for chromosome ${chr}."
      else
      echo "Error: vcf filtration for ${Output_file} failed"
    fi
endtime=`date +%s`
echo "chromosome ${chr} for ${input_cohort}  processed in $(($endtime - $starttime)) seconds."
      done
 
done
