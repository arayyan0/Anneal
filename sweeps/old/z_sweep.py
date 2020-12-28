#file: z_sweep.py
#author: Ahmed Rayyan
#date: June 22 2020 (updated)
#brief: creates a sweep varying Gz vs Kz
import sys
from numpy import linspace
import os

args = sys.argv

type, sub, l1, l2, Tf_pow, MS_pow, DS_pow, versions = [int(x) for x in args[1:9]]
cluster = input('WARNING: You better have picked the right cluster!\nInput one of RH/RE, and one of 1/2: ')

command = './sim %i %i %i %i %i %i'%(sub, l1, l2, Tf_pow, MS_pow, DS_pow)

plist = linspace(0,0.5,10+1)

path = 'out/z_job'
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

F = open('z_job.lst','w+')
for v in range(1, versions+1):
    if not os.path.exists(path+'/v_%i'%v):
        os.makedirs(path+'/v_%i'%v)
    for p in plist:
        F.write(command+' %.3f '%p + ' > ' + path+'/p_%.3f_.out\n'%p)
F.close()

F = open('z_job.sh','w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
F.write('#SBATCH --time=02:00:00'+'\n')
F.write('#SBATCH --job-name=z_job'+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load gcc'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 80 --joblog z_job.out <z_job.lst'+'\n')
F.close()
