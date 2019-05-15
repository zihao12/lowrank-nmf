args <- commandArgs(trailingOnly=TRUE)
dataname <- args[1]
K = 200
# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)

# These variables specify the names of the output files.
out.dir		 <- "../bigdata"
u.out.file  <- sprintf("%s_u_rsvd.csv", dataname)
v.out.file <- sprintf("%s_v_rsvd.csv", dataname)
d.out.file <- sprintf("%s_d_rsvd.csv", dataname)

# SET UP ENVIRONMENT
# ------------------
# Load packages and function definitions.
library(readr)
library(rsvd)
source(file.path("..","code","misc.R"))

# LOAD GTEX DATA
# --------------
cat("Loading GTEx data.\n")
read.counts.file <- file.path(data.dir,read.counts.file)
counts <- read.csv.matrix(read.counts.file)
cat(sprintf("Loaded %d x %d count matrix.\n",nrow(counts),ncol(counts)))

# LOAD INITIAL ESTIMATES
# ----------------------

# Get the number of factors ("topics").
#K <- ncol(F0)

# RUN NMF OPTIMIZATION METHOD
# ---------------------------
cat("Fitting Poisson topic model using nnmf.\n")
timing <- system.time(
	fit <- rsvd(counts, K)
)

cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

cat("Compute loglikelihood\n")



# WRITE NNMF RESULTS TO FILE
# --------------------------
cat("Writing results to file.\n")
u.out.file  <- file.path(out.dir,u.out.file)
v.out.file  <- file.path(out.dir,v.out.file)
d.out.file  <- file.path(out.dir,d.out.file)
write_csv(as.data.frame(fit$u),u.out.file,col_names = FALSE)
write_csv(as.data.frame(fit$v),v.out.file,col_names = FALSE)
write_csv(as.data.frame(fit$d),d.out.file,col_names = FALSE)

# SESSION INFO
# ------------
cat("Session info:\n")
sessionInfo()
