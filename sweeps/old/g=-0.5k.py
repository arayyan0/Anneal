#file: g=-0.5k.py
#author: Ahmed Rayyan
#date: March 4, 2020
#brief: finds the order at g = -0.5k

import sys
from numpy import linspace
import os

args = sys.argv

order, versions = args[1], int(args[2])
T_f_power, MS_power, DS_power = int(args[3]),int(args[4]),int(args[5])

orders = {
    "4": [2, 2, 1],
    "6" : [2, 3, 1],
    "12": [4, 3, 1],
    "18": [2, 3, 3],
    "30": [2, 5, 3],
    "50": [2, 5, 5]
}

sublattice, l1, l2 = orders[order]

p = 0.148
g = 0.5
a = 0.0

command = './sim %i %i %i %i %i %i'%(sublattice, l1, l2, T_f_power, MS_power, DS_power)

print("order = %s-site"%order)

path = "out/o_%s"%(order)
if not os.path.exists(path):
    os.makedirs(path)

F = open('o_%s.lst'%(order),'w+')
for version in range(1,versions+1):
    F.write(command+' %.3f %.3f %.3f'%(p,g,a) + ' > ' + path+'/v_%i.out\n'%(version))
F.close()
#
F = open('o_%s.sh'%(order),'w+')
F.write('#!/bin/bash\n\n')
F.write('#SBATCH --nodes=1\n')
F.write('#SBATCH --cpus-per-task=80\n')
F.write('#SBATCH --time=02:20:00\n')
F.write('#SBATCH --job-name=o_%s.lst\n\n'%(order))
F.write('cd $SLURM_SUBMIT_DIR\n\n')
F.write('module purge\n')
F.write('module load gnu-parallel\n')
F.write('parallel --jobs 0 --joblog o_%s.out < o_%s.lst\n'%(order,order))
F.close()
