#!/bin/bash

# common files, tools and parameters to all imputation scripts
WD=/path/to/Project
Bcftools=bcftools
Tabix=tabix
minimac4=minimac4
cohorts="${WD}/Cohorts_Info.csv"
fasta=/path/to/human_g1k_v37.fasta
eagle=/path/to/Eagle_v2.4.1/eagle
gmap=/path/to/genetic_map_hg19_withX.txt.gz
cpus=12
rounds=5
INFO="${WD}/extract_info_VCF.pl"
VCFtools=vcftools





