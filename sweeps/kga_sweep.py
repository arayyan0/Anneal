from sweep_lib import KGammaAnisotropyJobs

# min, max, spacing
p_val_list = [0.02, 0.02, 0.01, 'p']
a_val_list = [0, 1, 0.02, 'a']

# hc_or_kek, type, s, l1, l2
cluster_list = [0, 2, 6, 1, 1]

# Tf, MS_pow, DS_pow
param_list = [200, 4, 4]

run = 2
versions = 1

KGammaAnisotropyJobs('a', p_val_list, a_val_list, cluster_list, param_list, run, versions)
