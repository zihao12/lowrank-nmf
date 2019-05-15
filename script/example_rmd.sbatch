#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --exclusive
#SBATCH --time=16:00:00

# For reproducibility, I add version numbers to the "module load" calls.
module load pandoc
module load R/3.5.1

#Rscript -e "require ('knitr'); knit ('test.Rmd')"
Rscript -e "library(workflowr);wflow_build('../analysis/prepare_gtex_sample.Rmd',local = TRUE,view = FALSE,verbose = TRUE)"

