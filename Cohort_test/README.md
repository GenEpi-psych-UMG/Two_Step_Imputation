# Test dataset
Cohort test data for testing the provided workflow script for two-step imputation.
The folder includes VCF files of chromosome 10 and 11 of non-overlapping samples devided into 3 groups, genotype data was obtained from 1000G project (http://csg.sph.umich.edu/abecasis/MaCH/download/). The path to the test data directories in Cohorts_Info.csv needs be modified as specified in the workflow readme file.

Approximate runtime for the whole workflow using the test data is 20-30 minutes, using a 2.8 GHz Xeon CPU. The number of parallel CPU cores used for phasing and imputation can be increased to reduce the runtime. The number if used CPUs is set to 12 CPUs as specified in the configuration file.
