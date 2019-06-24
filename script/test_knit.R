args = commandArgs(trailingOnly=TRUE)

dataname  = args[1]
K = as.integer(args[2])

print(dataname)
print(K)
config = list(dataname = dataname, K = K)
# out.dir = "../docs"
out.dir = "../analysis"
rmd.dir = "../analysis"

rmdfile = "Investigate_pvalue_betanmf_log.Rmd"
rmdfile = file.path(rmd.dir, rmdfile)
outfile = sprintf("Investigate_pvalue_betanmf_%s_log_K_%d.html",config$dataname, config$K)
outfile = file.path(out.dir, outfile)


library(rmarkdown)
rmarkdown::render(rmdfile, params = config, output_file = outfile)
