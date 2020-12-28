#file: line_path.py
#author: Ahmed Rayyan
#date: August 20 2019
#brief: creating a line cutting through phase diagram

import numpy as np

vi = np.array([0.012, 0])
vf = np.array([0.025, 0.0500])

N = 100
ds = 1/N
n = 1
slist = np.linspace(0,1,N+1)

print("s = \n", slist)

command = './syed.out 1 0 0 0 0 0'
path = 'out' #make sure to create out/ folder first

F = open('ds-%.3f-%.0f.lst'%(ds,n),'w+')
for s in slist:
    p = (1-s)*vi[0] + s * vf[0]
    h = (1-s)*vi[1] + s * vf[1]
    F.write(command+' '+'%.5f'%p+' '+'%.5f'%h+' > '+path+'/s-%.3f.out'%s+'\n')
F.close()

#write batch file.

F = open('ds-%.3f-%.0f.sh'%(ds,n),'w+')
F.write('#!/bin/bash'+'\n'+'\n')
F.write('#SBATCH --nodes=1'+'\n')
F.write('#SBATCH --cpus-per-task=80'+'\n')
F.write('#SBATCH --time=3:13:10'+'\n')
F.write('#SBATCH --job-name=ds-%.3f'%ds+'\n'+'\n')
F.write('cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
F.write('module purge'+'\n')
F.write('module load intel/2018.2'+'\n')
F.write('module load gnu-parallel'+'\n')
F.write('parallel --jobs 1 --joblog ds-%.3f-%.0f.out < ds-%.3f-%.0f.lst'%(ds,n,ds,n))
F.close()
