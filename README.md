# Two-Step imputation workflow

This pipleline is designed to generate the intermediate panel that is required as an input for the two-step imputation workflow. The main workflow aims to overcome the technical bias due two batch effect between different genotyped cohorts. The whole worflow is represented in the figure below.
![Alt Text](Images/workflow.png)

Full description, testing and evaluation of the workflow is available in the manuscript.

## Table of Contents

- [Prerequisites](#Prerequisites)
- [Installation](#installation)
- [Outcome](#outcome)

### Prerequisites
- this project requires the follwoing softwares to be installed 
  [bcftools v1.14](https://samtools.github.io/bcftools/), [tabix 0.2.5](https://anaconda.org/bioconda/tabix), [minimac v4.1.4](https://genome.sph.umich.edu/wiki/Minimac4), [Eagle_v2.4.1](https://alkesgroup.broadinstitute.org/Eagle/), [VCFtools (0.1.17)](https://vcftools.sourceforge.net/index.html), [perl](https://www.perl.org/get.html)
- It also requires the presence of fasta file of human g1k (human_g1k_v37.fasta) and genetic_map_hg19 "genetic_map_hg19_withX.txt.gz"
Beside installation, the pathway to their location on the device needs to be declared in the configuration file as will be clarified later.

It also requires the genotyping data of the autosomal chromosomes for included cohorts in the two-step imputation to be one VCF file per chromomse, distinguished by numbers from 1-22. Each cohorts contains the genotyping data in a separate directory. For each cohort, the prefix and suffix of all chromosomes should be the same for each cohort.


## Installation
Clone the repository on your working directory, either by direct downlaod or using the git clone command
```bash
git clone https://github.com/username/repository.git
```

the repository includes stepwise scripts for the imputation workflow which require no further edits, only two files needs to be edited before starting the imputation:

1.  Cohorts_Info.csv
The table includes information about the genotyping cohorts, each row represents one cohorts including the pathway to the genotype data. example here for cohort "Cohort_A" where its VCF files in the following format:

```bash
bath/to/directory/Cohort_A_chr{i}.vcf.gz
```
each column should be filled for this cohort as follows;
Cohort = Cohort_A
Pathway= bath/to/directory
Prefix= Cohort_A_chr
Suffix= .vcf.gz

Note that the Pathway doesn't end with "/"

Save the table always incsv format

2. config.sh
Contains pathway and parameters needed for this pipeline, it is required to update the pathways of the working directory of the repository and other softwares in this file. there are other parameters incuded like number of used cpus and ietration rounds for imputation that can work without modifications but still modifiable according to the imputation's specific needs. 

## Usage

You can run the whole pipleine in one command 

which will take the working directory as an argument from the configuration source and run all commands of the pipeline at one time, after updating the cohort info and the configurations file just run the command:

```bash
bash bath/to/directory/Run_all.sh <bath/to/directory/config.sh>
```
as an alternative you can run each script separately as follows;

1. Alignment.sh
To insure the alignment of all included genotype data to the same strand, taking the location of fastafile as an input

```bash
bash bath/to/directory/Alignment.sh <bath/to/directory/config.sh>
```

2. vcf to bcf
Converting and indexing VCF files to BCF format for the follwoing phasing step.

```bash
bash bath/to/directory/VCF_BCF.sh <bath/to/directory/config.sh>
```

3. Phasing
phasing of the included genotype data
```bash
bash bath/to/directory/phase.sh <bath/to/directory/config.sh>
```

4. convert_refP
before the first round of imputation, a msav format of the phased genotyped data needs to be prepared using the follwoing command:


```bash
bash bath/to/directory/convert_RefP.sh <bath/to/directory/config.sh>
```
5.minimac4.sh
Running the imputation for all included cohorts, each cohort is imputed against all other included cohorts, all the output imputation files are saved in the same pathway of the imputed cohorts

```bash
bash bath/to/directory/minimac4.sh <bath/to/directory/config.sh>
```

6. R_Input.sh
The script generates an info file for each imputed genotpe file (one info file for each chromosme per cohort), containing variant information including R2 imputation quality. This should by used by the following R-scripted algorithm for building the intermediate panel.

```bash
bash bath/to/directory/R_Input.sh <bath/to/directory/config.sh>
```

7. R_Command.sh
This file runs the R algorithm (SNP_Selection.R) over the generated imputed gentoype files in step 5, the output of the algorithm is one filterout variants, which filter out variants from the corresponding imputed outcomes, and keep list files, which includes the variants to be filtered out from other imputed files. 

```bash
bash bath/to/directory/R_Command.sh <bath/to/directory/config.sh>
```

8. VCF_Filter_Vcftools.sh
This script uses VCFtools software to apply the filtration lists generated from the R-script to filter imputed genotype variants generated in step 5. 

```bash
bash bath/to/directory/VCF_Filter_Vcftools.sh <bath/to/directory/config.sh>
```

9. Merge_VCF_BCFtools.sh
The last step is to merge all the filtered files for each cohort into one vcf file per chromosome, the merging process is performed using BCF tools

```bash
bash bath/to/directory/Merge_VCF_BCFtools.sh <bath/to/directory/config.sh>
```


## Outcome
For each included cohort, an output folder is generated in the genotype data directory containg the intermediate panel files for autosomal chromosmes, ready to be used as an input for the second imputation using global reference panel (e.g. Michigan imputation server))
