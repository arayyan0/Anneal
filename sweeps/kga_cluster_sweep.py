from sweep_lib import SweepClusterJobs

p, a = 0.08, 0.00
cluster_info_list = [[2, 4, 1, n] for n in range(12,0,-1)]

# Tf, MS_pow, DS_pow
param_list = [200, 6, 4]

#ensure that num_anneal evenly divides cpus_per_task
cpus_per_task = 80
num_anneal = 20

run = 1

SweepClusterJobs(p, a, cluster_info_list, param_list, cpus_per_task, num_anneal, run)
