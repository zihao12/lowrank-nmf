## submit fitting job for algorithms in r
# "python fit_submit_r.py method data"

import os
import sys

data = sys.argv[1]
K = sys.argv[2]
K = int(K)

path = "./"

bshname = "fit_betanmf_{}_log_K{}.sbatch".format(data, K)
echo = "run fit_betanmf_log.R on {} with K {}".format(data,K)

output = '../output/fit_betanmf_{}_log_K{}.Rout'.format(data, K)
command = 'Rscript fit_betanmf_log.R {} {} > {}'.format(data,K, output)

with open(path + bshname, 'w') as rsh:
	with open(path + "example_R.sbatch", "r") as exa:
		for item in exa.readlines():
			rsh.write(item)
	rsh.write("echo '{}' \n".format(echo))
	rsh.write(command)

## submit job
print("sbatch {}".format(bshname))
os.system("sbatch {}".format(bshname))
