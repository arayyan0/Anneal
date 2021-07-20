from sweep_lib import SweepTemperatureJobs


#temperature array
T_min = 0.05
T_max = 1.5

#specifying step size
#dT = 0.05

#specifying number of steps
Nsteps = 12
dT = (T_max-T_min)/(Nsteps-1)

#T_min, #T_max, dT, label, spacing ('ari' or 'geo')
temp_list = [T_min, T_max, dT, 'Temp', 'geo_leftdense']


#l1, l2
cluster_list = [36, 36]

#ising_y, defect, lengthscale, number of defects
ham_list = [-0.0, -2.0, 1.5, 9]

run = 3
versions = 1

SweepTemperatureJobs(temp_list, cluster_list, ham_list, run, versions)
