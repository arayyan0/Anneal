#file: kg_sweep.py
#author: Ahmed Rayyan
#date: December 3 2019
#brief: explain after the fact
import sys
import numpy as np
import os
#
args = np.array(sys.argv)

sublattice, l1, l2 = [int(x) for x in args[1:4]]
final_T_power, max_msweeps_power, max_dsweeps_power = [int(x) for x in args[4:7]]
g, p = [float(x) for x in args[7:9]]
versions = int(args[9])

#alist = np.linspace(-0.1, 1, 22+1)
#alist = np.linspace(-0.1, 0.5, 12+1)
alist = np.linspace(-0.5, 1.0, 15+1)
print("a:")
print(alist)

command = "./sim %i %i %i %i %i %i"%(sublattice, l1, l2, final_T_power, max_msweeps_power, max_dsweeps_power)
# sunlattice, l1, l2, final_T_power, max_msweeps_power, max_dsweeps_power

path1 = "out/g_%.3f_p_%.3f/"%(g,p)
if not os.path.exists(path1):
    os.makedirs(path1)

F = open(path1+"SA_param.lst","w+")
F.write("type: 0\n")
F.write("(s, l1, l2) = (%i, %i, %i)\n"%(sublattice, l1, l2))
F.write("Tf = %.14f\n"%(0.9**final_T_power))
F.write("Metropolis sweeps = %i\n"%(10**max_msweeps_power))
F.write("Deterministic sweeps = %i\n"%(10**max_dsweeps_power))
F.close()

F = open("g_%.3f_p_%.3f.lst"%(g,p),"w+")
for version in range(1,versions+1):
    path = path1 + "v_%i"%(version)
    if not os.path.exists(path):
        os.makedirs(path)
    for a in alist:
        F.write(command+' %.3f %.3f %.3f'%(p,g,a) + ' > ' + path+'/a_%.3f_.out\n'%a)
F.close()

F = open("g_%.3f_p_%.3f.sh"%(g,p),'w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
F.write('#SBATCH --time=04:30:00'+'\n')
F.write('#SBATCH --job-name=g_%.3f_p_%.3f'%(g,p)+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 80 --joblog g_%.3f_p_%.3f.out < g_%.3f_p_%.3f.lst\n'%(g,p,g,p))
F.close()


# sublattice, l1, l2 = [int(x) for x in args[1:4]]
# final_T_power, max_msweeps_power, max_dsweeps_power = [int(x) for x in args[4:7]]
# g, a = [float(x) for x in args[7:9]]
# versions = int(args[9])
#
# plist = np.linspace(0, 0.5, 50+1)
# print("p:")
# print(plist)
#
# command = "./sim %i %i %i %i %i %i"%(sublattice, l1, l2, final_T_power, max_msweeps_power, max_dsweeps_power)
# # sunlattice, l1, l2, final_T_power, max_msweeps_power, max_dsweeps_power
#
# path1 = "out/g_%.3f_a_%.3f/"%(g,a)
# if not os.path.exists(path1):
#     os.makedirs(path1)
#
# F = open(path1+"SA_param.lst","w+")
# F.write("type: 0\n")
# F.write("(s, l1, l2) = (%i, %i, %i)\n"%(sublattice, l1, l2))
# F.write("Tf = %.14f\n"%(0.9**final_T_power))
# F.write("Metropolis sweeps = %i\n"%(10**max_msweeps_power))
# F.write("Deterministic sweeps = %i\n"%(10**max_dsweeps_power))
# F.close()
#
# F = open("g_%.3f_a_%.3f.lst"%(g,a),"w+")
# for version in range(1,versions+1):
#     path = path1 + "v_%i"%(version)
#     if not os.path.exists(path):
#         os.makedirs(path)
#     for p in plist:
#         F.write(command+' %.3f %.3f %.3f'%(p,g,a) + ' > ' + path+'/p_%.3f_.out\n'%(p))
# F.close()
#
# F = open("g_%.3f_a_%.3f.sh"%(g,a),'w+')
# F.write('#!/bin/bash'+'\n'+'\n')
# F.write('#SBATCH --nodes=1'+'\n')
# F.write('#SBATCH --cpus-per-task=80'+'\n')
# F.write('#SBATCH --time=01:00:00'+'\n')
# F.write('#SBATCH --job-name=g_%.3f_a_%.3f'%(g,a)+'\n'+'\n')
# F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
# F.write('module purge'+'\n')
# F.write('module load gnu-parallel'+'\n')
# F.write('parallel --jobs 80 --joblog g_%.3f_a_%.3f.out < g_%.3f_a_%.3f.lst\n'%(g,a,g,a))
# F.close()
