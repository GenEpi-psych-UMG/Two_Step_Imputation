#!/bin/bash

# Source the configuration file to load paths and parameters
source "$1"

# Check if the necessary variables are loaded properly
if [ -z "$Bcftools" ] || [ -z "$VCFtools" ] || [ -z "$cohorts" ]; then
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

  for j in "${!cohort_names[@]}"; do
    target_cohort="${cohort_names[$j]}"
    target_pathway="${cohort_paths[$j]}"

# Trim any leading/trailing whitespace from variables
  input_pathway=$(echo "$input_pathway" | sed 's/[[:space:]]*$//')
  target_pathway=$(echo "$target_pathway" | sed 's/[[:space:]]*$//')
  input_cohort=$(echo "$input_cohort" | xargs)
  target_cohort=$(echo "$target_cohort" | xargs)

    # Ensure input cohort and target cohort are not the same
    if [ "$input_cohort" != "$target_cohort" ]; then

      # Loop through chromosomes 1 to 2
      for chr in {1..22}; do

        # Define the input, target, and output file paths
        Input_file="${input_pathway}/${input_cohort}_chr${chr}_imputed_${target_cohort}.vcf.gz"
        filter_file="${input_pathway}/Keep_list_main_chr${chr}_${target_cohort}.txt"
        Output_file="${input_pathway}/${input_cohort}_chr${chr}_imputed_${target_cohort}"
        WD="${input_pathway}"

        # Run the minimac4 command
echo "Input_file= '$Input_file'"
echo "filter_file= '$filter_file'"
echo "Output_file= '$Output_file'"

 # Change to the working directory
        if cd "${WD}"; then
          echo "Changed to working directory: ${WD}"
        else
          echo "Error: Could not change to working directory: ${WD}"
          exit 1
        fi

# Check if the filter file exists in the working directory
if [ -r "${filter_file}" ]; then
    # If the file exists, run the first vcftools command
    "$VCFtools" --gzvcf "${Input_file}" \
         --chr ${chr} \
         --positions "${filter_file}" \
         --recode \
         --out "${Output_file}_addon"
bgzip -c "${Output_file}_addon.recode.vcf" > "${Output_file}_addon.recode.vcf.gz"
"$Tabix" -p vcf "${Output_file}_addon.recode.vcf.gz"
else
    # If the file does not exist, run the alternative vcftools command
        "$VCFtools" --gzvcf "${Input_file}" \
         --chr ${chr} \
         --exclude-positions "${target_pathway}/filterout_main_chr${chr}.txt" \
         --recode \
         --out "${Output_file}_kept"
bgzip -c "${Output_file}_kept.recode.vcf" > "${Output_file}_kept.recode.vcf.gz"
"$Tabix" -p vcf "${Output_file}_kept.recode.vcf.gz"

fi

if [ $? -eq 0 ]; then
      echo "vcf filtration for of ${Output_file} is done"
    else
      echo "Error: vcf filtration for ${Output_file} failed"
    fi
      done
    fi
  done
done
