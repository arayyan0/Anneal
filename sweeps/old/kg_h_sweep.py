#file: kg_h_sweep.py
#author: Ahmed Rayyan
#date: June 18 2020 (updated)
#brief: creates a field sweep at a given point in the KG phase space
import sys
from numpy import linspace
import os

args = sys.argv

type, sub, l1, l2, Tf_pow, MS_pow, DS_pow, versions, hth = [int(x) for x in args[1:10]]
p, gp = [float(x) for x in args[10:12]]
cluster = input('WARNING: You better have picked the right cluster!\nInput one of RH/RE, and one of 1/2: ')

command = './sim %i %i %i %i %i %i'%(sub, l1, l2, Tf_pow, MS_pow, DS_pow)

a = float(args[12])
hlist = linspace(0, 1.5, 30+1)
print("h:")
print(hlist)

path = 'out/hth_%i_p_%.3f_gp_%.3f/a_%.3f'%(hth,p,gp,a)
if not os.path.exists(path):
  os.makedirs(path)
#
F = open(path+"/SA_param.lst",'w+')
F.write("type: %i\n"%type)
F.write("(s, l1, l2) = (%i, %i, %i)\n"%(sub, l1, l2))
F.write("Tf = %.20f\n"%(0.9)**Tf_pow)
F.write("Metropolis sweeps = %i\n"%(10)**MS_pow)
F.write("Deterministic sweeps = %i\n"%(10)**DS_pow)
F.write("Cluster: %s\n"%cluster)
F.close()

F = open('a_%.3f.lst'%a,'w+')
for v in range(1, versions+1):
    if not os.path.exists(path+'/v_%i'%v):
        os.makedirs(path+'/v_%i'%v)
    for h in hlist:
        F.write(command+' %i %.3f %.3f %.3f %.3f'%(hth,p,gp,a,h) + ' > ' + path+'/v_%i/a_%.3f_h_%.3f_.out\n'%(v,a,h))
F.close()

F = open('a_%.3f.sh'%a,'w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
F.write('#SBATCH --time=02:00:00'+'\n')
F.write('#SBATCH --job-name=a_%.3f'%a+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load gcc'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 80 --joblog a_%.3f.out < a_%.3f.lst'%(a,a)+'\n')
F.close()
