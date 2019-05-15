## submit fitting job for matlab algorithms
# "python fit_submit_matlab.py method data"

import os
import sys

method = sys.argv[1]
data = sys.argv[2]

path = "./"

bshname = "fit_{}_{}.sbatch".format(data, method)
echo = "run fit_{}.m on {} data".format(method, data)

matcom1 = "dataname = '{}'".format(data)
matcom2 = "run('fit_{}.m')".format(method)
matcomm = ' "{};{}; exit;"'.format(matcom1, matcom2)

output = '../output/fit_{}_{}.out'.format(data, method)

command = 'matlab -nodisplay -nosplash -nodesktop -r {} > {}'.format(matcomm, output)

with open(path + bshname, 'w') as rsh:
	with open(path + "example_matlab.sbatch", "r") as exa:
		for item in exa.readlines():
			rsh.write(item)
	rsh.write("echo '{}' \n".format(echo))
	rsh.write(command)

## submit job
print("sbatch {}".format(bshname))
os.system("sbatch {}".format(bshname))
