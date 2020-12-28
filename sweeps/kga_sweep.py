from sweep_lib import KGammaAnisotropyJobs

# min, max, spacing
p_val_list.append('p') = [0.25, 0.50, 0.01, "p"]
a_val_list.append('a') = [0.00, 0.00, 1, "a"]

# hc_or_kek, type, s, l1, l2
cluster_list = [0, 2, 2, 3, 3]

# Tf, MS_pow, DS_pow
param_list = [200, 3, 4]

run = 1

KGammaAnisotropyJobs('p', p_val_list, a_val_list, cluster_list, param_list, run)
