import itertools as it
import numpy as np
import os

class SweepPhaseDiagramJobs:
    def __init__(self, cluster_list, params_list, params_label_list, sa_list,
                 cpus_per_task, num_anneal, run, versions, node):

        self.ClusterList = cluster_list

        self.Tf_Pow, self.MS_Pow, self.DS_Pow = sa_list

        self.JobTitle = f"jobrun_{run}"
        print(self.JobTitle)
        self.OutputPath = f"out/{self.JobTitle}"
        if not os.path.exists(self.OutputPath):
            os.makedirs(self.OutputPath)

        self.Labels = params_label_list
        self.Params = []
        for params in params_list:
            x_min, x_max, dx = params
            if abs(x_max - x_min) < pow(10,-8):
                self.Params.append(np.array([x_min]))
            else:
                N_x = round(1+ (x_max-x_min)/dx)
                self.Params.append(np.linspace(x_min, x_max, N_x))
        print(self.ClusterList)
        for i in range(len(self.Labels)):
            print(f"{self.Labels[i]}:\n{self.Params[i]}")

        self.WriteJobDescription()
        print("Writing .sh file...")
        print("Please ensure that you changed the number of nodes and time from the default!")
        self.WriteSHFile(cpus_per_task,num_anneal,node)
        self.WriteLSTFile(num_anneal, versions)

    def WriteJobDescription(self):
        F = open(self.OutputPath+"/SA_param.lst", 'w+')
        for i in range(len(self.Labels)):
            F.write(f"{self.Labels[i]}-array: {self.Params[i]}\n")
        F.write("all files in this folder have the following global parameters.\n")
        F.write(f"Tf = {(0.9)**self.Tf_Pow:.20f}\n")
        F.write(f"Metropolis sweeps = {(10)**self.MS_Pow}\n")
        F.write(f"Deterministic sweeps = {(10)**self.DS_Pow}\n")
        F.close()

    def WriteSHFile(self,cpus_per_task,num_anneal,node):
        F = open(f'{self.JobTitle}.sh','w+')
        F.write(f'#!/bin/bash'+'\n'+'\n')
        F.write(f'#SBATCH --nodes={node}'+'\n')
        F.write(f'#SBATCH --ntasks-per-node={cpus_per_task}'+'\n')
        F.write(f'#SBATCH --time=00:30:00'+'\n')
        F.write(f'#SBATCH --job-name={self.JobTitle}'+'\n'+'\n')
        F.write(f'cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
        F.write(f'module purge'+'\n')
        F.write(f'module load gcc openmpi gnu-parallel'+'\n')
        #F.write(f'parallel --jobs {int(node*cpus_per_task/num_anneal)} --joblog {self.OutputPath}/{self.JobTitle}.out'+
        #          f' < {self.JobTitle}.lst\n')
        F.write(f'parallel --joblog {self.OutputPath}/{self.JobTitle}.out --wd $PWD'+
                  f' < {self.JobTitle}.lst\n')
        F.close()

    def WriteLSTFile(self, num_anneal, versions):
        # for lat, shape, l1, l2, l3 in self.ClusterList:
        File = open(f'{self.JobTitle}.lst','w+')
        commandbegin = f"mpirun -np {num_anneal} --bind-to none ./sim "
        commandend = f"{self.Tf_Pow} {self.MS_Pow} {self.DS_Pow}"
        a = [list(param) for param in self.Params]
        product =   list(it.product(*a))
        for v in range(1,versions+1):
            for lat, shape, l1, l2, l3 in self.ClusterList:
                if not os.path.exists(self.OutputPath+f'/lat_{lat}_s_{shape}_l1_{l1}_l2_{l2}_l3_{l3}/v_{v}'):
                    os.makedirs(self.OutputPath+f'/lat_{lat}_s_{shape}_l1_{l1}_l2_{l2}_l3_{l3}/v_{v}/')
                for p in product:
                    param_text_list = ''.join(map(str, [f'{p[i]:.6f} ' for i in range(len(self.Params))]))
                    param_label_list = ''.join(map(str, [f'{self.Labels[i]}_{p[i]:.6f}_' for i in range(len(self.Params))]))
                    File.write(commandbegin + \
                               f"{lat} {shape} {l1} {l2} {l3} " + \
                               param_text_list + \
                               commandend + ' > ' + f'{self.OutputPath}'+ f'/lat_{lat}_s_{shape}_l1_{l1}_l2_{l2}_l3_{l3}/v_{v}' + \
                               '/'+ param_label_list+'.out\n')
        File.close()
