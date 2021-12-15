from sweep1_lib import SweepPhaseDiagramJobs
import numpy as np


# lattice, shape, l1, l2, l3
cluster_list = [
                [1,1,3,1,1],
                [1,0,1,6,1],
                [1,0,1,4,1],
                [1,0,1,3,1],
                [1,2,1,1,1],
                [1,0,1,2,1],
                [1,0,1,1,1],
               ]
# IF USING MULTIPLE CLUSTERS, ENSURE THE LARGER ONES GET DONE FIRST

# parameter entry: min, max, spacing
entry = 0
if entry == 0:
    #jtau2 + jq2 + jo2 = 1
    t_val_list,  t_label  = [ 0.0075, 0.150, 0.0075], "t"
    p_val_list,  p_label  = [ 0.490, 0.490, 0.005], "p"
    params_list = [t_val_list, p_val_list]
    params_label_list = [t_label, p_label]

elif entry == 1:
    #abs(jtau) =1
    jq_val_list, jq_label = [0.000, 0.000, 0.001], "jq"
    jo_val_list, jo_label = [0.000, 0.600, 0.060], "jo"
    params_list = [jq_val_list, jo_val_list]
    params_label_list = [jq_label, jo_label]

jb_val_list, jb_label = [0.000, 0.000, 0.001], "jb"
# jb_val_list, jb_label = [np.sqrt(2.0/5.0), np.sqrt(2.0/5.0), 0.001], "jb"
h_val_list,  h_label  = [0.000, 0.000, 0.001], "h"

params_list = params_list + [jb_val_list, h_val_list]
params_label_list = params_label_list + [jb_label,h_label]


# need to implement defect_list
# defect_list = [0, 0, 1]

# Tf, MS_pow, DS_pow
sa_list = [200, 5, 4]

#ensure that num_anneal divides cpus_per_task
cpus_per_task = 80
num_anneal = 4

run = 1
versions = 1

SweepPhaseDiagramJobs(cluster_list, params_list, params_label_list, sa_list,
                      cpus_per_task, num_anneal, run, versions)
