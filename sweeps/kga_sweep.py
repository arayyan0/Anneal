from sweep_lib import KGammaAnisotropyJobs

# min, max, spacing
p_val_list = [0.08, 0.08, 0.01, 'p']
a_val_list = [0.00, 0.00, 0.01, 'a']

# hc_or_kek, type, s, l1, l2
cluster_list = [0, 2, 4, 1, 1]

# Tf, MS_pow, DS_pow
param_list = [200, 5, 4]

run = 1
versions = 1

KGammaAnisotropyJobs('p', p_val_list, a_val_list, cluster_list, param_list, run, versions)
