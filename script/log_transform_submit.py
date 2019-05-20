## submit fitting job for algorithms in r
# "python fit_submit_r.py method data"

import os
import sys

data = sys.argv[1]

path = "./"

bshname = "log_transform_{}.sbatch".format(data)
echo = "run log_transform.R on {} data".format(data)

output = '../output/log_transform_{}.Rout'.format(data)
command = 'Rscript log_transform.R {} > {}'.format(data, output)

with open(path + bshname, 'w') as rsh:
	with open(path + "example_R.sbatch", "r") as exa:
		for item in exa.readlines():
			rsh.write(item)
	rsh.write("echo '{}' \n".format(echo))
	rsh.write(command)

## submit job
print("sbatch {}".format(bshname))
os.system("sbatch {}".format(bshname))
