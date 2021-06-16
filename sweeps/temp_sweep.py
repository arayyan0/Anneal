from sweep_lib import SweepTemperatureJobs

#temperature array
temp_list = [2, 6, 0.25, 'Temp']

#l1, l2
cluster_list = [18, 18]

#ising_y, defect
ham_list = [-1, 0]

run = 3
versions = 5

SweepTemperatureJobs(temp_list, cluster_list, ham_list, run, versions)
