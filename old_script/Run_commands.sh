#!/bin/bash

# Usage: bash run_pipeline.sh /path/to/config.sh

# Check if the config file is provided
if [ -z "$1" ]; then
  echo "Error: Path to config.sh is required."
  exit 1
fi

# Source the configuration file
CONFIG_FILE="$1"
source "$CONFIG_FILE"

# Check if WD is set in the config file
if [ -z "$WD" ]; then
  echo "Error: 'WD' variable not set in config.sh."
  exit 1
fi

# Step 1: Alignment
echo "Running Alignment..."
bash "$WD/Alignment.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: Alignment.sh failed."
  exit 1
fi

# Step 2: VCF to BCF Conversion and Indexing
echo "Converting VCF to BCF..."
bash "$WD/VCF_BCF.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: VCF_BCF.sh failed."
  exit 1
fi

# Step 3: Phasing
echo "Running Phasing..."
bash "$WD/phase.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: phase.sh failed."
  exit 1
fi

# Step 4: Convert to msav format (Reference Panel preparation)
echo "Converting to RefP format..."
bash "$WD/convert_RefP.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: convert_RefP.sh failed."
  exit 1
fi

# Step 5: Imputation using Minimac4
echo "Running Minimac4 Imputation..."
bash "$WD/minimac4.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: minimac4.sh failed."
  exit 1
fi

# Step 6: Generate R input files
echo "Generating R Input..."
bash "$WD/R_Input.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: R_Input.sh failed."
  exit 1
fi

# Step 7: Run the R Algorithm for SNP Selection
echo "Running R Command..."
bash "$WD/R_Command.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: R_Command.sh failed."
  exit 1
fi

# Step 8: Apply VCF Filtering using VCFtools
echo "Applying VCF Filtering..."
bash "$WD/VCF_Filter_Vcftools.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: VCF_Filter_Vcftools.sh failed."
  exit 1
fi

# Step 9: Merge filtered VCFs using BCFtools
echo "Merging VCF Files..."
bash "$WD/Merge_VCF_BCFtools.sh" "$CONFIG_FILE"
if [ $? -ne 0 ]; then
  echo "Error: Merge_VCF_BCFtools.sh failed."
  exit 1
fi

echo "Pipeline completed successfully."
