args <- commandArgs(trailingOnly=TRUE)
dataname <- args[1]

# method = "rough"
# dataname = "test"

# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)


out.dir		 <- "../output"
hist.out.file  <- sprintf("%s_hist.png", dataname)

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

## turn count to log (add e)
e = 10^(-10)
counts = max(counts, e)


cat("plot and save.\n")
# These variables specify the names of the output files.


hist.out.file  <- file.path(out.dir,hist.out.file)
png(hist.out.file)
hist(log10(counts), breaks = 100, xlab = "log of counts", main = sprintf("counts (log) of %s", dataname))


# SESSION INFO
# ------------
cat("Session info:\n")
sessionInfo()
