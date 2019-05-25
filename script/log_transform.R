args <- commandArgs(trailingOnly=TRUE)
dataname <- args[1]
# SCRIPT SETTINGS
# ---------------
# These variables specify the names of the input files.
data.dir           <- "../bigdata"
read.counts.file   <- sprintf("%s.csv", dataname)


# These variables specify the names of the output files.
out.dir		 <- "../bigdata"
logcounts.out.file  <- sprintf("%s_log.csv", dataname)

# SET UP ENVIRONMENT
# ------------------
# Load packages and function definitions.
library(readr)
source(file.path("..","code","misc.R"))

# LOAD GTEX DATA
# --------------
cat("Loading count data.\n")
read.counts.file <- file.path(data.dir,read.counts.file)
counts <- read.csv.matrix(read.counts.file)
cat(sprintf("Loaded %d x %d count matrix.\n",nrow(counts),ncol(counts)))


# DO LOG TRANSFORM
counts_log = logtrans(counts)


# WRITE NNMF RESULTS TO FILE
# --------------------------
cat("Writing results to file.\n")
logcounts.out.file  <- file.path(out.dir,logcounts.out.file)
write_csv(as.data.frame(counts_log),logcounts.out.file,col_names = FALSE)

# SESSION INFO
# ------------
cat("Session info:\n")
sessionInfo()
