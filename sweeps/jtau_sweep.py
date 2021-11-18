from sweep1_lib import SweepPhaseDiagramJobs
import numpy as np

# lattice, shape, l1, l2, l3
cluster_list = [1, 0, 6, 6, 1]

# min, max, spacing
entry = 0
if entry == 0:
    #jtau2 + jq2 + jo2 +jb2 = 1
    p1_val_list, p1_label = [0, 0.5, 0.025],    "t1"
    p2_val_list, p2_label = [0.5, 0.5, 0],      "t2"
    p3_val_list, p3_label = [0.75, 0.75, 0.01], "p"
elif entry == 1:
    #jtau =1, p1 = jb_unitless, p2 = jq_unitless, p3 = jo_unitless
    p1_val_list, p1_label = [0, 0.0, 0.025],    "jb"
    p2_val_list, p2_label = [0.0, 0.0, 0],      "jq"
    p3_val_list, p3_label = [0.0, 0.0, 0.01],   "jo"

h_val_list,      h_label      = [0.0, 0.0, 0], "h"

# need to implement defect_list
# defect_list = [0, 0, 1]

params_list = np.array([p1_val_list, p2_val_list, p3_val_list, h_val_list])
params_label_list = [p1_label, p2_label, p3_label, h_label]

# Tf, MS_pow, DS_pow
sa_list = [50, 3, 3]

#ensure that num_anneal divides cpus_per_task
cpus_per_task = 6
num_anneal = 6

run = 1

SweepPhaseDiagramJobs(cluster_list, params_list, params_label_list, sa_list,
                      cpus_per_task, num_anneal, run)
