#file: kg_a_sweep.py
#author: Ahmed Rayyan
#date: January 28 2020
#brief: creates anisotropic +G->-K sweep

import numpy as np
import sys
from numpy import linspace
import os

args = sys.argv
type, sub, l1, l2, Tf_pow, MS_pow, DS_pow, versions = [int(x) for x in args[1:9]]
cluster = input('WARNING: You better have picked the right cluster!\nInput one of RH/RE, and one of 1/2: ')
run=int(args[10])
command = './sim %i %i %i %i %i %i'%(sub, l1, l2, Tf_pow, MS_pow, DS_pow)

# g = 0.5
# a = float(args[9])
# print("a:")
# print(a)
#
# plst = linspace(0, 0.5, 25+1)
# print("p:")
# print(plst)
#
# path = 'out/g_%.3f_a_%.3f'%(g,a)
# if not os.path.exists(path):
#   os.makedirs(path)
#
# F = open(path+"/SA_param.lst",'w+')
# F.write("type: %i\n"%type)
# F.write("(s, l1, l2) = (%i, %i, %i)\n"%(sub, l1, l2))
# F.write("Tf = %.20f\n"%(0.9)**Tf_pow)
# F.write("Metropolis sweeps = %i\n"%(10)**MS_pow)
# F.write("Deterministic sweeps = %i\n"%(10)**DS_pow)
# F.write("Cluster: %s\n"%cluster)
# F.close()
#
#
# F = open('g_%.3f_a_%.3f.lst'%(g,a),'w+')
# for v in range(1, versions+1):
#     if not os.path.exists(path+'/v_%i'%v):
#         os.makedirs(path+'/v_%i'%v)
#     for p in plst:
#         F.write(command+' %.3f %.3f %.3f'%(p,g,a) + ' > ' + path+'/v_%i/p_%.3f_.out\n'%(v,p))
# F.close()
# #
# F = open('g_%.3f_a_%.3f.sh'%(g,a),'w+')
# F.write('#!/bin/bash'+'\n'+'\n')
# F.write('#SBATCH --nodes=1'+'\n')
# F.write('#SBATCH --cpus-per-task=80'+'\n')
# F.write('#SBATCH --time=02:00:00'+'\n')
# F.write('#SBATCH --job-name=g_%.3f_a_%.3f_'%(g,a)+'\n'+'\n')
# F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
# F.write('module purge'+'\n')
# F.write('module load gcc'+'\n')
# F.write('module load gnu-parallel'+'\n')
# F.write('parallel --jobs 80 --joblog g_%.3f_a_%.3f.out < g_%.3f_a_%.3f.lst'%(g,a,g,a)+'\n')
# F.close()

g = 0.5
p = float(args[9])
print("p:")
print(p)
# alst = linspace(-0.5, 1, 30+1)
# #####
# alst = linspace(0.4,0.65,25+1)
# #####
# alst = linspace(-0.03,0.21,80+1) #80+1 or 24+1
# alst = linspace(0.1,0.5, ) #80+1 or 25+1
# alst = linspace(0.4,0.8, ) #80+1 or 25+1
# alst = linspace(0,0.4,20+1)
# alst = np.array([0.00,0.10,0.45,0.60])
# alst = np.array([0.999])
alst = [0.8]
print("a:")
print(alst)

path = 'out/g_%.3f_p_%.3f_r_%i'%(g,p,run)
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

F = open('g_%.3f_p_%.3f_r_%i.lst'%(g,p,run),'w+')
for v in range(1, versions+1):
    if not os.path.exists(path+'/v_%i'%v):
      os.makedirs(path+'/v_%i'%v)
    for a in alst:
        F.write(command+' %.3f %.3f %.3f'%(p,g,a) + ' > ' + path+'/v_%i/a_%.3f_.out\n'%(v,a))
F.close()
#
F = open('g_%.3f_p_%.3f_r_%i.sh'%(g,p,run),'w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
F.write('#SBATCH --time=00:30:00'+'\n')
F.write('#SBATCH --job-name=g_%.3f_p_%.3f_r_%i'%(g,p,run)+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load gcc'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 80 --joblog out/g_%.3f_p_%.3f_r_%i/g_%.3f_p_%.3f_r_%i.out < g_%.3f_p_%.3f_r_%i.lst'%(g,p,run,g,p,run,g,p,run)+'\n')
F.close()
