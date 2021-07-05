from sweep_lib import SweepTemperatureJobs


#temperature array
T_min = 0.05
T_max = 0.55

#specifying step size
#dT = 0.05

#specifying number of steps
Nsteps = 12
dT = (T_max-T_min)/(Nsteps-1)

#T_min, #T_max, dT, label, spacing ('ari' or 'geo')
temp_list = [T_min, T_max, dT, 'Temp', 'geo']


#l1, l2
cluster_list = [36, 36]

#ising_y, defect, number of defects
ham_list = [-0.4, -2.0, 9]

run = 1
versions = 1

SweepTemperatureJobs(temp_list, cluster_list, ham_list, run, versions)
