from sweep_lib import SweepTemperatureJobs

#temperature array
#T_min, #T_max, dT, label, spacing ('ari' or 'geo')
temp_list = [0.05, 0.6, 0.05, 'Temp', 'geo']

#l1, l2
cluster_list = [36, 36]

#ising_y, defect
ham_list = [-0.0, -0.0]

run = 1
versions = 3

SweepTemperatureJobs(temp_list, cluster_list, ham_list, run, versions)
