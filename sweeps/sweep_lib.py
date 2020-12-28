import numpy as np
import os
import itertools as it

class SweepPhaseDiagramJobs:
    '''
    '''
    def __init__(self, x_val_list, y_val_list, cluster_list, param_list, run):
        '''
        Parameters
        x_val_list (list): x-values parameters
        y_val_list (list): y-values parameters. if only interested in 1d at some
                           fixed y = y_fixed, use [y_fixed, y_fixed, 1].
        s, l1, l2   (int):
        '''
        self.XArray, self.YArray = self.DiscretizedGrid(x_val_list[:3], y_val_list[:3])
        self.XLabel, self.YLabel = x_val_list[3], y_val_list[3]

        print(f"{self.XLabel}:\n{self.XArray}")
        print(f"{self.YLabel}:\n{self.YArray}")

        self.HcOrKek, self.ClusterType, self.S, self.L1, self.L2 = cluster_list
        self.Tf_Pow, self.MS_Pow, self.DS_Pow = param_list

        self.JobTitle = f"jobrun_{run}"
        self.OutputPath = f"out/{self.JobTitle}"
        if not os.path.exists(self.OutputPath):
            os.makedirs(self.OutputPath)

        self.WriteJobDescription()
        print("Writing .sh file...")
        print("Please ensure that you changed the number of nodes and time from the default!")
        self.WriteSHFile()

    def WriteJobDescription(self):
        F = open(self.OutputPath+"/SA_param.lst", 'w+')
        F.write(f"{self.XLabel}-array: {self.XArray}\n")
        F.write(f"{self.YLabel}-array: {self.YArray}\n")
        F.write("all files in this folder have the following global parameters.\n")
        F.write(f"hc or kekule: {self.HcOrKek}\n")
        F.write(f"type: {self.ClusterType}\n")
        F.write(f"(s, l1, l2) = ({self.S}, {self.L1}, {self.L2})\n")
        F.write(f"Tf = {(0.9)**self.Tf_Pow:.20f}\n")
        F.write(f"Metropolis sweeps = {(10)**self.MS_Pow}\n")
        F.write(f"Deterministic sweeps = {(10)**self.DS_Pow}\n")
        F.close()

    def WriteSHFile(self):
        F = open(f'{self.JobTitle}.sh','w+')
        F.write(f'#!/bin/bash'+'\n'+'\n')
        F.write(f'#SBATCH --nodes=1'+'\n')
        F.write(f'#SBATCH --cpus-per-task=80'+'\n')
        F.write(f'#SBATCH --time=00:30:00'+'\n')
        F.write(f'#SBATCH --job-name={self.JobTitle}'+'\n'+'\n')
        F.write(f'cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
        F.write(f'module purge'+'\n')
        F.write(f'module load gcc'+'\n')
        F.write(f'module load gnu-parallel'+'\n')
        F.write(f'parallel --jobs 80 --joblog {self.OutputPath}/{self.JobTitle}.out'+
                  f' < {self.JobTitle}.lst\n')
        F.close()

    def DiscretizedGrid(self, x_val_list, y_val_list):
        '''
        Parameters
        x_val_list (list-like): x values to be discretized
        y_val_list (list-like): y values to be discretized
        Returns
        x_array (numpy.ndarray): discretized x-values
        y_array (numpy.ndarray): discretized y-values
        '''
        x_min, x_max, dx = x_val_list
        y_min, y_max, dy = y_val_list
        N_x  = int(1 + (x_max-x_min)/dx)
        N_y  = int(1 + (y_max-y_min)/dy)
        x_array = np.linspace(x_min, x_max, N_x)
        y_array = np.linspace(y_min, y_max, N_y)
        return x_array, y_array

class KGammaAnisotropyJobs(SweepPhaseDiagramJobs):
    def __init__(self, which_swept, p_list, a_list, cluster_list, param_list, run, versions):
        if which_swept == 'p':
            super().__init__(p_list, a_list, cluster_list, param_list, run)
            self.pArray, self.aArray = self.XArray, self.YArray
        elif which_swept == 'a':
            super().__init__(a_list, p_list, cluster_list, param_list, run)
            self.aArray, self.pArray = self.XArray, self.YArray
        self.WriteLSTFile(versions)

    def WriteLSTFile(self, versions):
        command = f"./sim {self.ClusterType} {self.S} {self.L1} {self.L2} " +\
                     f"{self.Tf_Pow} {self.MS_Pow} {self.DS_Pow}"

        product = list(it.product(list(self.pArray), list(self.aArray)))
        File = open(f'{self.JobTitle}.lst','w+')
        for v in range(1, versions+1):
            if not os.path.exists(self.OutputPath+f'/v_{v}'):
                os.makedirs(self.OutputPath+f'/v_{v}')
            for prod in product:
                File.write(command+f' {prod[0]:.3f} 0.5 {prod[1]:.3f}' + ' > ' +
                        self.OutputPath+f'/v_{v}/p_{prod[0]:.3f}_a_{prod[1]:.3f}_.out\n')
        File.close()
