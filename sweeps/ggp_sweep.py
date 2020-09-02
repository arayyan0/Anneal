#file: kg_a_sweep.py
#author: Ahmed Rayyan
#date: January 28 2020
#brief: creates anisotropic +G->-K sweep

import sys
from numpy import linspace
import os

args = sys.argv
type, sub, l1, l2, Tf_pow, MS_pow, DS_pow, versions = [int(x) for x in args[1:9]]
cluster = input('WARNING: You better have picked the right cluster!\nInput one of RH/RE, and one of 1/2: ')

command = './sim %i %i %i %i %i %i'%(sub, l1, l2, Tf_pow, MS_pow, DS_pow)

plst = linspace(0.4, 0.6, 40+1)
print("p:")
print(plst)

path = 'out/ggp'
if not os.path.exists(path):
  os.makedirs(path)

F = open(path+"/SA_param.lst",'w+')
F.write("type: %i\n"%type)
F.write("(s, l1, l2) = (%i, %i, %i)\n"%(sub, l1, l2))
F.write("Tf = %.20f\n"%(0.9)**Tf_pow)
F.write("Metropolis sweeps = %i\n"%(10)**MS_pow)
F.write("Deterministic sweeps = %i\n"%(10)**DS_pow)
F.write("Cluster: %s\n"%cluster)
F.close()

F = open('ggp.lst','w+')
for v in range(1, versions+1):
    if not os.path.exists(path+'/v_%i'%v):
      os.makedirs(path+'/v_%i'%v)
    for p in plst:
        F.write(command+' %.3f'%(p) + ' > ' + path+'/v_%i/p_%.3f_.out\n'%(v,p))
F.close()
#
F = open('ggp.sh','w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
F.write('#SBATCH --time=03:20:00'+'\n')
F.write('#SBATCH --job-name=ggp'+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load gcc'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 80 --joblog ggp.out < ggp.lst'+'\n')
F.close()
