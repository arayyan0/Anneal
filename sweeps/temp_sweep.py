from sweep_lib import SweepTemperatureJobs


#temperature array
T_min = 0.05
T_max = 0.6

#specifying step size
#dT = 0.05

#specifying number of steps
Nsteps = 80
dT = (T_max-T_min)/(Nsteps-1) 

#T_min, #T_max, dT, label, spacing ('ari' or 'geo')
temp_list = [T_min, T_max, dT, 'Temp', 'geo']


#l1, l2
cluster_list = [36, 36]

#ising_y, defect
ham_list = [-0.0, -0.0]

run = 1
versions = 3

SweepTemperatureJobs(temp_list, cluster_list, ham_list, run, versions)
