args <- commandArgs(trailingOnly=TRUE)
method <- args[1]
dataname <- args[2]

# method = "rough"
# dataname = "test"

# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)

f.file  <- sprintf("%s_factors_%s.csv", dataname, method)
l.file  <- sprintf("%s_loadings_%s.csv", dataname, method)





# SET UP ENVIRONMENT
# ------------------
# Load packages and function definitions.
library(readr)
source(file.path("..","code","misc.R"))
# source(file.path("..","code","utility.R"))
# source(file.path("..","code","betanmf.R"))

# LOAD GTEX DATA
# --------------
cat("Loading GTEx data.\n")
read.counts.file <- file.path(data.dir,read.counts.file)
counts <- read.csv.matrix(read.counts.file)
cat(sprintf("Loaded %d x %d count matrix.\n",nrow(counts),ncol(counts)))

# LOAD  ESTIMATES
# ----------------------
cat("Loading approximation.\n")
f.file  <- file.path(data.dir,f.file)
l.file  <- file.path(data.dir,l.file)


F                 <- read.csv.matrix(f.file)
L                 <- read.csv.matrix(l.file)
K = ncol(F)

Lam <- L %*% t(F)


cat("plot and save.\n")
# These variables specify the names of the output files.
out.dir		 <- "../output"
resid.out.file  <- sprintf("%s_%s_K%d_resid_pearson.png", dataname, method,K)
pval.out.file <- sprintf("%s_%s_K%d_pval_pearson.png",dataname, method, K)
pval2.out.file <- sprintf("%s_%s_K%d_pval_poisson.png",dataname, method, K)
resid.out.file  <- file.path(out.dir,resid.out.file)
pval.out.file  <- file.path(out.dir,pval.out.file)
pval2.out.file  <- file.path(out.dir,pval2.out.file)



## compute and plot resid and pval 
main = sprintf("%s on %s with K = %d", method, dataname, K)
pval_resid_plot_pearson(counts, Lam, resid.out.file, pval.out.file,  main)
pval_plot_poisson(counts, Lam, pval2.out.file,  main)

# SESSION INFO
# ------------
cat("Session info:\n")
sessionInfo()
