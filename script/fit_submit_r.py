## submit fitting job for algorithms in r
# "python fit_submit_r.py method data"

import os
import sys

method = sys.argv[1]
data = sys.argv[2]

path = "./"

bshname = "fit_{}_{}.sbatch".format(data, method)
echo = "run fit_{}.R on {} data".format(method, data)

output = '../output/fit_{}_{}.Rout'.format(data, method)
command = 'Rscript fit_{}.R {} > {}'.format(method, data, output)

with open(path + bshname, 'w') as rsh:
	with open(path + "example_R.sbatch", "r") as exa:
		for item in exa.readlines():
			rsh.write(item)
	rsh.write("echo '{}' \n".format(echo))
	rsh.write(command)

## submit job
print("sbatch {}".format(bshname))
os.system("sbatch {}".format(bshname))
