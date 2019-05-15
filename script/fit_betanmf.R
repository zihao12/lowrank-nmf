args <- commandArgs(trailingOnly=TRUE)
dataname <- args[1]
# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)
init.factors.file  <- sprintf("%s_factors_rough.csv", dataname)
init.loadings.file <- sprintf("%s_loadings_rough.csv", dataname)

cat(sprintf("%s\n",read.counts.file))
cat(sprintf("%s\n",init.factors.file))

# These variables specify the names of the output files.
out.dir		 <- "../bigdata"
factors.out.file  <- sprintf("%s_factors_betanmfr.csv", dataname)
loadings.out.file <- sprintf("%s_loadings_betanmfr.csv", dataname)

# SET UP ENVIRONMENT
# ------------------
# Load packages and function definitions.
library(readr)
library(NNLM)
source(file.path("..","code","misc.R"))
source(file.path("..","code","utility.R"))
source(file.path("..","code","betanmf.R"))

# LOAD GTEX DATA
# --------------
cat("Loading GTEx data.\n")
read.counts.file <- file.path(data.dir,read.counts.file)
counts <- read.csv.matrix(read.counts.file)
cat(sprintf("Loaded %d x %d count matrix.\n",nrow(counts),ncol(counts)))

# LOAD INITIAL ESTIMATES
# ----------------------
cat("Loading initial estimates of factors and loadings.\n")
init.factors.file  <- file.path(data.dir,init.factors.file)
init.loadings.file <- file.path(data.dir,init.loadings.file)
F0                 <- read.csv.matrix(init.factors.file)
L0                 <- read.csv.matrix(init.loadings.file)
cat(sprintf("Loaded %d x %d factors matrix, ",nrow(F0),ncol(F0)))
cat(sprintf("and %d x %d loadings matrix.\n",nrow(L0),ncol(L0)))

# Get the number of factors ("topics").
K <- ncol(F0)

# RUN NMF OPTIMIZATION METHOD
# ---------------------------
cat("Fitting Poisson topic model using nnmf.\n")
timing <- system.time(
	fit <- betanmf(counts, L0, t(F0), numiter=20)
)

cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

cat("Compute loglikelihood\n")
F <- t(fit$B)
L <- fit$A

out <- compute_ll(t(counts),F,t(L))
cat(sprintf("method type: %s\n
        poisson_ll :%0.12f\n
        multinom_ll:%0.12f\n",out$type, out$pois_ll,out$multinom_ll))


# WRITE NNMF RESULTS TO FILE
# --------------------------
cat("Writing results to file.\n")
factors.out.file  <- file.path(out.dir,factors.out.file)
loadings.out.file <- file.path(out.dir,loadings.out.file)
write_csv(as.data.frame(F),factors.out.file,col_names = FALSE)
write_csv(as.data.frame(L),loadings.out.file,col_names = FALSE)

# SESSION INFO
# ------------
cat("Session info:\n")
sessionInfo()
