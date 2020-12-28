#file: kg_a_sweep.py
#author: Ahmed Rayyan
#date: January 28 2020
#brief: creates anisotropic +G->-K sweep

import sys
from numpy import linspace
import os

args = sys.argv

order = args[1]
T_f_power, MS_power, DS_power = int(args[2]),int(args[3]),int(args[4])

orders = {
    "4": [2, 2, 1],
    "6" : [2, 3, 1],
    "12": [4, 3, 1],
    "18": [2, 3, 3],
    "30": [2, 5, 3],
    "50": [2, 5, 5]
}

sublattice, l1, l2 = orders[order]

command = './sim %i %i %i %i %i %i'%(sublattice, l1, l2, T_f_power, MS_power, DS_power)
#l1, l2, sublattice, T_f_power, max_msweeps_power, max_dsweeps_power
g = 0.5
alist = linspace(0, 0.4, 10+1)
plist = linspace(0, 0.26, 13+1)

print("order = %s-site"%order)
print("p:")
print(plist)
print("g = %.3f"%g)
print("a:")
print(alist)

for version in range(1,10+1):
    path = "out/o_%s/v_%s"%(order,version)
    if not os.path.exists(path):
        os.makedirs(path)

    F = open('o_%s_v_%i.lst'%(order,version),'w+')
    for p in plist:
      for a in alist:
        F.write(command+' %.3f %.3f %.3f'%(p,g,a) + ' > ' + path+'/p_%.3f_a_%.3f_.out\n'%(p, a))
    F.close()

    F = open('o_%s_v_%i.sh'%(order,version),'w+')
    F.write('#!/bin/bash\n\n')
    F.write('#SBATCH --nodes=1\n')
    F.write('#SBATCH --cpus-per-task=80\n')
    F.write('#SBATCH --time=02:20:00\n')
    F.write('#SBATCH --job-name=o_%s_v_%i.lst\n\n'%(order,version))
    F.write('cd $SLURM_SUBMIT_DIR\n\n')
    F.write('module purge\n')
    F.write('module load gnu-parallel\n')
    F.write('parallel --jobs 0 --joblog o_%s_v_%i.out < o_%s_v_%i.lst\n'%(order, version, order, version))
    F.close()
