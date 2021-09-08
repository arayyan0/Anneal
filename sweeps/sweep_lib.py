import numpy as np
import os
import itertools as it

class SweepClusterJobs:
    '''
    '''
    def __init__(self, p, a, cluster_info_list, param_list, cpus_per_task, num_anneal, run):
        self.JobTitle = f"jobrun_{run}"
        print(self.JobTitle)
        self.OutputPath = f"out/{self.JobTitle}"
        if not os.path.exists(self.OutputPath):
            os.makedirs(self.OutputPath)


        print("(p,a): " + f"({p:.3f},{a:.3f})")
        print(cluster_info_list)

        self.WriteJobDescription(p, a, cluster_info_list, param_list)
        print("Writing .sh file...")
        print("Please ensure that you changed the number of nodes and time from the default!")
        self.WriteSHFile(cpus_per_task, num_anneal)
        self.WriteLSTFile(p, a, cluster_info_list, param_list, num_anneal)


    def WriteJobDescription(self, p, a, cluster_info_list, param_list):
        F = open(self.OutputPath+"/SA_param.lst", 'w+')
        F.write("cluster list:\n")
        for cluster in cluster_info_list:
            F.write(f"[{cluster[0]}, {cluster[1]}, {cluster[2]}, {cluster[3]}]\n")
        F.write("all files in this folder have the following global parameters.\n")
        F.write(f"p, a: {p}, {a}\n")
        F.write(f"Tf = {(0.9)**param_list[0]:.20f}\n")
        F.write(f"Metropolis sweeps = {(10)**param_list[1]}\n")
        F.write(f"Deterministic sweeps = {(10)**param_list[2]}\n")
        F.close()

    def WriteSHFile(self,cpus_per_task,num_anneal):
        F = open(f'{self.JobTitle}.sh','w+')
        F.write(f'#!/bin/bash'+'\n'+'\n')
        F.write(f'#SBATCH --nodes=1'+'\n')
        F.write(f'#SBATCH --ntasks-per-node={cpus_per_task}'+'\n')
        F.write(f'#SBATCH --time=00:30:00'+'\n')
        F.write(f'#SBATCH --job-name={self.JobTitle}'+'\n'+'\n')
        F.write(f'cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
        F.write(f'module purge'+'\n')
        F.write(f'module load gcc openmpi gnu-parallel'+'\n')
        F.write(f'parallel --jobs {int(cpus_per_task/num_anneal)} --joblog {self.OutputPath}/{self.JobTitle}.out'+
                  f' < {self.JobTitle}.lst\n')
        F.close()

    def WriteLSTFile(self, p, a, cluster_info_list, param_list, num_anneal):
        commandbegin = f"mpirun -np {num_anneal} --bind-to none ./sim "
        commandend = f" {param_list[0]} {param_list[1]} {param_list[2]}"
        File = open(f'{self.JobTitle}.lst','w+')
        for cluster in cluster_info_list:
            File.write(commandbegin + \
                       f'{cluster[0]} {cluster[1]} {cluster[2]} {cluster[3]} {p} 0.5 {a}' + \
                       commandend + \
                       ' > ' + \
                       f'{self.OutputPath}'+ \
                       f'/ct_{cluster[0]}_s_{cluster[1]}_l1_{cluster[2]}_l2_{cluster[3]}_.out\n')
        File.close()


class SweepPhaseDiagramJobs:
    '''
    '''
    def __init__(self, x_val_list, y_val_list, cluster_list, param_list,
                 cpus_per_task, num_anneal, run):
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
        print(self.JobTitle)
        self.OutputPath = f"out/{self.JobTitle}"
        if not os.path.exists(self.OutputPath):
            os.makedirs(self.OutputPath)

        self.WriteJobDescription()
        print("Writing .sh file...")
        print("Please ensure that you changed the number of nodes and time from the default!")
        self.WriteSHFile(cpus_per_task,num_anneal)

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

    def WriteSHFile(self,cpus_per_task,num_anneal):
        F = open(f'{self.JobTitle}.sh','w+')
        F.write(f'#!/bin/bash'+'\n'+'\n')
        F.write(f'#SBATCH --nodes=1'+'\n')
        F.write(f'#SBATCH --ntasks-per-node={cpus_per_task}'+'\n')
        F.write(f'#SBATCH --time=00:30:00'+'\n')
        F.write(f'#SBATCH --job-name={self.JobTitle}'+'\n'+'\n')
        F.write(f'cd $SLURM_SUBMIT_DIR'+'\n'+'\n')
        F.write(f'module purge'+'\n')
        F.write(f'module load gcc openmpi gnu-parallel'+'\n')
        F.write(f'parallel --jobs {int(cpus_per_task/num_anneal)} --joblog {self.OutputPath}/{self.JobTitle}.out'+
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
        N_x  = round(1 + (x_max-x_min)/dx)
        N_y  = round(1 + (y_max-y_min)/dy)
        x_array = np.linspace(x_min, x_max, N_x)
        y_array = np.linspace(y_min, y_max, N_y)
        return x_array, y_array

class KGammaAnisotropyJobs(SweepPhaseDiagramJobs):
    def __init__(self, which_swept, p_list, a_list, cluster_list, param_list,
                 cpus_per_task, num_anneal, run):
        if which_swept == 'p':
            super().__init__(p_list, a_list, cluster_list, param_list,
                             cpus_per_task, num_anneal, run)
            self.pArray, self.aArray = self.XArray, self.YArray
        elif which_swept == 'a':
            super().__init__(a_list, p_list, cluster_list, param_list,
                             cpus_per_task, num_anneal, run)
            self.aArray, self.pArray = self.XArray, self.YArray
        self.WriteLSTFile(num_anneal)

    def WriteLSTFile(self, num_anneal):
        commandbegin = f"mpirun -np {num_anneal} --bind-to none ./sim "
        commandend = f" {self.Tf_Pow} {self.MS_Pow} {self.DS_Pow}"
        product = list(it.product(list(self.pArray), list(self.aArray)))
        File = open(f'{self.JobTitle}.lst','w+')
        for prod in product:
            File.write(commandbegin + \
                       f'{self.ClusterType} {self.S} {self.L1} {self.L2} {prod[0]:.3f} 0.5 {prod[1]:.3f}' + \
                       commandend + \
                       ' > ' + \
                       f'{self.OutputPath}'+ \
                       f'/p_{prod[0]:.3f}_a_{prod[1]:.3f}_.out\n')
        File.close()

class SweepTemperatureJobs:
    def __init__(self, temp_list, cluster_list, ham_list, run, versions):
        x_min, x_max, dx = temp_list[:3]
        N_x = round(1+ (x_max-x_min)/dx)
        self.XLabel = temp_list[3]
        spacing = temp_list[4]

        if spacing == 'ari':
            self.XArray = np.linspace(x_min, x_max, N_x)
        elif spacing == 'geo_rightdense':
            self.XArray = (x_max+x_min) - np.geomspace(x_min, x_max, N_x)
        elif spacing == 'geo_leftdense':
	        self.XArray = np.geomspace(x_min, x_max, N_x)

        print(f"{self.XLabel}:\n{self.XArray}")

        self.L1, self.L2 = cluster_list

        self.JQuad, self.JOcto, self.DefectQuad, self.DefectOcto, self.Lengthscale, self.NumDefects = ham_list

        self.JobTitle = f"jobrun_{run}"
        print(self.JobTitle)
        self.OutputPath = f"out/{self.JobTitle}"
        if not os.path.exists(self.OutputPath):
            os.makedirs(self.OutputPath)

        self.WriteJobDescription()
        print("Writing .sh file...")
        print("Please ensure that you changed the number of nodes and time from the default!")
        self.WriteSHFile()

        self.WriteLSTFile(versions)

    def WriteJobDescription(self):
        F = open(self.OutputPath+"/global_param.lst", 'w+')
        F.write(f"{self.XLabel}-array: {self.XArray}\n")
        F.write("all files in this folder have the following global parameters.\n")
        F.write(f"(l1, l2) = ({self.L1}, {self.L2})\n")
        F.write(f"j_quad, j_octo, defect_q, defect_o, lengthscale, num_defects = {self.JQuad} {self.JOcto} {self.DefectQuad} {self.DefectOcto}, {self.Lengthscale}, {self.NumDefects}\n")
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
        F.write(f'module load gcc openmpi gnu-parallel'+'\n')
        F.write(f'parallel --jobs 80 --joblog {self.OutputPath}/{self.JobTitle}.out'+
                  f' < {self.JobTitle}.lst\n')
        F.close()

    def WriteLSTFile(self, versions):
        command = f"./sim {self.L1} {self.L2} {self.NumDefects} {self.JQuad} {self.JOcto} {self.DefectQuad} {self.DefectOcto} {self.Lengthscale}"

        File = open(f'{self.JobTitle}.lst','w+')
        for v in range(1, versions+1):
            if not os.path.exists(self.OutputPath+f'/v_{v}'):
                os.makedirs(self.OutputPath+f'/v_{v}')
            for x in self.XArray:
                File.write(command+f' {x:.6f}' + ' > ' +
                        self.OutputPath+f'/v_{v}/{self.XLabel}_{x:.6f}_.out\n')
        File.close()
