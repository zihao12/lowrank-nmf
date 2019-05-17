args <- commandArgs(trailingOnly=TRUE)
dataname <- args[1]
method <- args[2]
d <- as.integer(args[3])
# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)

u.file  <- sprintf("%s_u_%s.csv", dataname, method)
v.file  <- sprintf("%s_v_%s.csv", dataname, method)
d.file  <- sprintf("%s_d_%s.csv", dataname, method)



# These variables specify the names of the output files.
out.dir		 <- "../output"
hist.out.file  <- sprintf("%s_approx_%s_d%d_hist.png", dataname, method, d)
scatter.out.file <- sprintf("%s_approx_%s_d%d_scatter.png",dataname, method, d)

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
u.file  <- file.path(data.dir,u.file)
v.file  <- file.path(data.dir,v.file)
d.file  <- file.path(data.dir,d.file)


U                 <- read.csv.matrix(u.file)
V                 <- read.csv.matrix(v.file)
D                 <- read.csv.matrix(d.file)

countshat <- U[,1:d] %*% diag(D[1:d]) %*% t(V[,1:d])


cat(sprintf("Reconstruct %d x %d count matrix.\n",nrow(countshat),ncol(countshat)))

cat("plot hist of absolute difference.\n")
hist.out.file  <- file.path(out.dir,hist.out.file)
png(hist.out.file)
hist(abs(counts - countshat), main = "hist of abs(X - Xhat)")


cat("scatter plot of X vs Xhat.\n")
scatter.out.file  <- file.path(out.dir,scatter.out.file)
png(scatter.out.file)
plot(counts, countshat, main = "X vs Xhat", xlab = "counts", ylab = "countshat")

# SESSION INFO
# ------------
cat("Session info:\n")
sessionInfo()
