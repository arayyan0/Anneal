from sweep_lib import SweepTemperatureJobs


#temperature array
T_min = 0.05
T_max = 1.5

#specifying step size
#dT = 0.05

#specifying number of steps
Nsteps = 92
dT = (T_max-T_min)/(Nsteps-1)

#T_min, #T_max, dT, label, spacing ('ari' or 'geo')
temp_list = [T_min, T_max, dT, 'Temp', 'geo_leftdense']


#l1, l2
cluster_list = [36, 36]

#j_quad, j_octo, defect_q, defect_o, lengthscale, num_defects
ham_list = [-0.0,-0.0,-0.0,-0.0, 2.0, 1]

run = 1
versions = 2

SweepTemperatureJobs(temp_list, cluster_list, ham_list, run, versions)
