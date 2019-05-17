#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --exclusive
#SBATCH --time=16:00:00

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1

Rscript visual_approx.R gtex rsvd 50
