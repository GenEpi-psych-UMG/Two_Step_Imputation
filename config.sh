#!/bin/bash

# common files, tools and parameters to all imputation scripts
WD=/daten/gtnas/users/nasrk/Imputation/Scripts/Project
Bcftools=bcftools
Tabix=tabix
minimac4=minimac4
cohorts=/daten/gtnas/users/nasrk/Imputation/Scripts/Project/Cohorts_Info.csv
fasta=/daten/gtnas/users/nasrk/Imputation/Data/human_g1k_v37.fasta
eagle=/daten/gtnas/users/nasrk/Imputation/Programs/Eagle_v2.4.1/eagle
gmap=/daten/gtnas/users/nasrk/Imputation/Programs/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz
cpus=12
rounds=5
INFO=/daten/gtnas/users/nasrk/Imputation/Scripts/Project/extract_info_VCF.pl
VCFtools=vcftools





