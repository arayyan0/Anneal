from sweep1_lib import SweepPhaseDiagramJobs
import numpy as np

# lattice, shape, l1, l2, l3
cluster_list = [1, 0, 2, 2, 1]

# min, max, spacing
theta1_val_list, theta1_label = [0, 0.5, 0.025], "t1"
theta2_val_list, theta2_label = [0.5, 0.5, 0], "t2"
phi_val_list,    phi_label    = [0.75, 0.75, 0.01], "p"
h_val_list,      h_label      = [0.0, 0.0, 0], "h"

params_list = np.array([theta1_val_list, theta2_val_list, phi_val_list, h_val_list])
params_label_list = [theta1_label, theta2_label, phi_label, h_label]

# Tf, MS_pow, DS_pow
sa_list = [50, 3, 3]

#ensure that num_anneal divides cpus_per_task
cpus_per_task = 80
num_anneal = 16

run = 1

SweepPhaseDiagramJobs(cluster_list, params_list, params_label_list, sa_list,
                      cpus_per_task, num_anneal, run)
