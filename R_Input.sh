#!/bin/bash

# Source the configuration file to load paths and parameters
source "$1"

# Check if the necessary variables are loaded properly
if [ -z "$INFO" ] ||  [ -z "$cohorts" ]; then
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

# Create a separate output log file for each input cohort
  output_log="${input_pathway}/${input_cohort}_INFO_files.txt"
  echo -e "Path\tReference\tchr" > "$output_log"

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
        Output_file="${input_pathway}/${input_cohort}_chr${chr}_imputed_${target_cohort}"


echo "Output_file= '$Output_file'"

# Run the minimac4 command
"$INFO" ${input_pathway}/ ${input_cohort}_chr${chr}_imputed_${target_cohort}
if [ $? -eq 0 ]; then
      echo "imputation of ${Output_file} is done"
# Append the output file path and target cohort (reference) to the cohort-specific log file
          echo -e ""${Output_file}.all.info.txt"\t${target_cohort}\t${chr}" >> "$output_log"

    else
      echo "Error: imputation of ${Output_file} failed"
    fi
      done
    fi
  done
done
