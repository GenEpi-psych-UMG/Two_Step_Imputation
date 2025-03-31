#!/bin/bash

# Source the configuration file to load paths and parameters
source "$1"

# Check if the necessary variables are loaded properly
if [ -z "$minimac4" ] || [ -z "$cpus" ] || [ -z "$cohorts" ]; then
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
        Ref_file="${target_pathway}/${target_cohort}_chr${chr}_phased_RefP.msav"
        Input_file="${input_pathway}/${input_cohort}_chr${chr}_phased.vcf.gz"
        Output_file="${input_pathway}/${input_cohort}_chr${chr}_imputed_${target_cohort}"

        # Run the minimac4 command
echo "Input_file= '$Input_file'"
echo "Ref_file= '$Ref_file'"
echo "Output_file= '$Output_file'"

# Run the minimac4 command
        "$minimac4" ${Ref_file} \
                  ${Input_file} \
                  --output ${Output_file}.vcf.gz \
                  -f GT \
                  --all-typed-sites \
                  --threads ${cpus}
if [ $? -eq 0 ]; then
      echo "imputation of ${Output_file} is done"
    else
      echo "Error: imputation of ${Output_file} failed"
    fi
      done
    fi
  done
done
