#file: kkj_sweep.py
#author: Ahmed Rayyan
#date: January 23 2020
#brief: creates a sweep through the KKJ phase line, includes niagara job script
import sys
import numpy as np
import os

args = sys.argv
p_i, p_f, N, name = float(args[1]), float(args[2]), int(args[3]), args[4]
plist = np.linspace(p_i, p_f, N+1)
#plist = np.array([0.065,0.195,0.02,0.03])
#plist3 = np.array([0.05,0.06,0.07])

#for i in range(0,5):
#    nplist=(0.04*i)+plist3
#    plist = np.concatenate((plist, nplist))

print("p:")
print(plist)

command = './sim 12 80 4 5'
#length, T_f_power, max_msweeps_power, max_dsweeps_power

path = 'out'
if not os.path.exists(path):
    os.mkdir(path)

F = open('%s.lst'%name,'w+')
for p in plist:
    F.write(command+' '+'%.3f'%p + ' > ' + path+'/p_%.3f_.out\n'%p)
F.close()

F = open('%s.sh'%name,'w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
#F.write('#SBATCH --time='+time+'\n')
F.write('#SBATCH --time=01:00:00'+'\n')
F.write('#SBATCH --job-name=%s'%name+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load gcc'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 0 --joblog %s.out < %s.lst\n'%(name,name))
F.close()
