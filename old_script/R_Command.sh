#!/bin/bash

# Source the configuration file to load paths and parameters
source "$1"

# Check if the necessary variables are loaded properly
if [ -z "$WD" ] ||  [ -z "$cohorts" ]; then
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

# Trim any leading/trailing whitespace from variables
  input_pathway=$(echo "$input_pathway" | sed 's/[[:space:]]*$//')
  target_pathway=$(echo "$target_pathway" | sed 's/[[:space:]]*$//')
  input_cohort=$(echo "$input_cohort" | xargs)
  target_cohort=$(echo "$target_cohort" | xargs)

            # Define the input, target, and output file paths
        Input_file="${input_pathway}/${input_cohort}_INFO_files.txt"
        
cd ${input_pathway}
Rscript "$WD"/SNP_Selection.R "$Input_file"
 
if [ $? -eq 0 ]; then
      echo "filtration lists for ${input_cohort} is generated"

    else
      echo "Error: filtration lists for ${input_cohort} failed"
    fi
      done
