args <- commandArgs(trailingOnly=TRUE)
dataname <- args[1]
K <-  as.integer(args[2])

# dataname = "test"
# K = 5
# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)
#init.factors.file  <- sprintf("%s_factors_rough.csv", dataname)
#init.loadings.file <- sprintf("%s_loadings_rough.csv", dataname)

# These variables specify the names of the output files.
out.dir		 <- "../bigdata"
factors.out.file  <- sprintf("%s_log_factors_betanmf_K%d.csv", dataname,K)
loadings.out.file <- sprintf("%s_log_loadings_betanmf_K%d.csv", dataname,K)

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
cat("Loading and transforming data.\n")
read.counts.file <- file.path(data.dir,read.counts.file)
counts <- read.csv.matrix(read.counts.file)
counts_log <- logtrans_add1(counts)
cat(sprintf("Loaded and transformed %d x %d count matrix.\n",nrow(counts_log),ncol(counts_log)))

# LOAD INITIAL ESTIMATES
# ----------------------
cat("Initial estimates of factors and loadings using NNMF.\n")
set.seed(123)
init = nnmf(counts_log,K,method = "scd",loss = "mkl",rel.tol = 1e-8,
              n.threads = 0,max.iter = 10,inner.max.iter = 4,trace = 1,
              verbose = 2)

# RUN NMF OPTIMIZATION METHOD
# ---------------------------
cat("Fitting Poisson topic model.\n")
timing <- system.time(
	fit <- betanmf(counts, init$W, init$H, numiter=20, verbose = T, eval_every = 5)
)
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

#cat("Compute loglikelihood\n")
F <- t(fit$B)
L <- fit$A

#out <- compute_ll(t(counts),F,t(L))
# cat(sprintf("method type: %s\n
#         poisson_ll :%0.12f\n
#         multinom_ll:%0.12f\n",out$type, out$pois_ll,out$multinom_ll))


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
