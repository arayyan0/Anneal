#file: kg_a_sweep.py
#author: Ahmed Rayyan
#date: January 28 2020
#brief: creates anisotropic +G->-K sweep

import numpy as np
import sys
import os

args = sys.argv
hc_or_kek, type, sub, l1, l2, Tf_pow, MS_pow, DS_pow, versions = [int(x) for x in args[1:10]]
run=int(args[10])
command = f"./sim {type} {sub} {l1} {l2} {Tf_pow} {MS_pow} {DS_pow}"

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

p = float(args[11])
print(f"p: {p}")
g = 0.5
# alst = np.linspace(-0.4, 0, 80+1)
alst = np.linspace(0.005, 0.805, 80+1)
print("a:")
print(alst)

path = f'out/p_{p:.3f}_r_{run}'
if not os.path.exists(path):
  os.makedirs(path)

F = open(path+"/SA_param.lst",'w+')
F.write(f"hc or kekule: {hc_or_kek}\n")
F.write(f"type: {type}\n")
F.write(f"(s, l1, l2) = ({sub}, {l1}, {l2})\n")
F.write(f"Tf = {(0.9)**Tf_pow:.20f}\n")
F.write(f"Metropolis sweeps = {(10)**MS_pow}\n")
F.write(f"Deterministic sweeps = {(10)**DS_pow}\n")
F.close()

F = open(f'p_{p:.3f}_r_{run}.lst','w+')
for v in range(1, versions+1):
    if not os.path.exists(path+f'/v_{v}'):
      os.makedirs(path+f'/v_{v}')
    for a in alst:
        F.write(command+f' {p:.3f} {g:.3f} {a:.3f}' + ' > ' + path+f'/v_{v}/a_{a:.3f}_.out\n')
F.close()
#
F = open(f'p_{p:.3f}_r_{run}.sh','w+')
F.write(f'#!/bin/bash'+'\n'+'\n')
F.write(f'#SBATCH --nodes=1'+'\n')
F.write(f'#SBATCH --cpus-per-task=80'+'\n')
F.write(f'#SBATCH --time=00:30:00'+'\n')
F.write(f'#SBATCH --job-name=p_{p}_r_{run}'+'\n'+'\n')
F.write(f'cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write(f'module purge'+'\n')
F.write(f'module load gcc'+'\n')
F.write(f'module load gnu-parallel'+'\n')
F.write(f'parallel --jobs 80 --joblog out/p_{p:.3f}_r_{run}/p_{p:.3f}_r_{run}.out < p_{p:.3f}_r_{run}.lst'+'\n')
F.close()
