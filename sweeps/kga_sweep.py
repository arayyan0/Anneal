from sweep_lib import KGammaAnisotropyJobs

# min, max, spacing
p_val_list = [0.073, 0.091, 0.002, 'p']
a_val_list = [0.0, 0.0, 0.01, 'a']

# hc_or_kek, type, s, l1, l2
cluster_list = [0, 2, 2, 1, 2]

# Tf, MS_pow, DS_pow
param_list = [200, 5, 4]

#ensure that num_anneal evenly divides cpus_per_task
cpus_per_task = 80
num_anneal = 16

run = 1

KGammaAnisotropyJobs('p', p_val_list, a_val_list, cluster_list, param_list,
                     cpus_per_task, num_anneal, run)
