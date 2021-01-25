# -*- coding: utf-8 -*-

import numpy as np
import sys,os
import pickle
import matplotlib.pyplot as plt
import random
import time

# Import tools
import tools

Network_params = {
    "exc_exc_recurrent": 0.178,
    "exc_inh_recurrent": 0.233,
    "inh_inh_recurrent": -2.70,
    "inh_exc_recurrent": -2.01
}

# New seed
random.seed(int(time.time()))

# Type of processing: either ratio of conductances (0) or chemogenetic
# manipulations (1)
proc_type = 0

# Manipulation variable
label_var = 'g'

# File ID
IDs = ["Sim0","Sim1"]

# External input rates
ext_rates = [1.5,2.0] # (spikes/s)

# Automatic search of g_ex and g_in, or the chemo. variable range
if proc_type == 0:
    g_ex_range = []
    g_in_range = []

cont_var_range = []

# Load all files
for id in IDs:
    if proc_type == 0:
        # Limits
        limits = [5.,15.]
        g_ex_range_aux = []
        g_in_range_aux = []
    else:
        cont_var_range_aux = []
        # Limits
        # limits = [-75.,-69.]
        limits = [-53.,-51.5]

    dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'results/'+id))
    ldir = os.listdir(dir)
    ldir.sort()

    for filename in ldir:
        if 'AMPA' in filename:
            # g_ex
            if proc_type == 0:
                m1 = filename.find("gex_")+4
                m2 = filename.find("_gin_")
                g_ex = float(filename[m1:m2])

                # g_in
                m1 = filename.find("gin_")+4
                m2 = filename.find(".AMPA")
                g_in = float(filename[m1:m2])

                # Ext freq
                m1 = filename.find("rate_")+5
                m2 = filename.find("_gex")
                v0 = float(filename[m1:m2])

                g_ratio = -g_in*Network_params["inh_exc_recurrent"]/(g_ex*Network_params["exc_exc_recurrent"])
                if (g_ex == 1.0 or g_in == 1.0) and (v0 in ext_rates) and (g_ratio>= limits[0] and g_ratio <=limits[1]):
                    g_ex_range_aux.append(g_ex)
                    g_in_range_aux.append(g_in)

            else:
                # cont_var
                m1 = filename.find(cont_var+"_")+len(cont_var)+1
                m2 = filename.find(".AMPA")
                var = float(filename[m1:m2])

                # Ext freq
                m1 = filename.find("rate_")+5
                m2 = filename.find("_"+cont_var)
                v0 = float(filename[m1:m2])

                if v0 in ext_rates and (var>= limits[0] and var <=limits[1]):
                    cont_var_range_aux.append(var)


    # Select only unique values
    if proc_type == 0:
        g_ex_range_aux = np.array(g_ex_range_aux)
        g_in_range_aux = np.array(g_in_range_aux)

        g_ex_range_aux = np.unique(g_ex_range_aux)
        g_in_range_aux = np.unique(g_in_range_aux)

        # Add to g_ex_range and g_in_range
        random.shuffle(g_ex_range_aux)
        random.shuffle(g_in_range_aux)
        g_ex_range.append(g_ex_range_aux)
        g_in_range.append(g_in_range_aux)
    else:
        cont_var_range_aux = np.array(cont_var_range_aux)
        cont_var_range_aux = np.unique(cont_var_range_aux)

        # Add to cont_var_range
        random.shuffle(cont_var_range_aux)
        cont_var_range.append(cont_var_range_aux)

# For conductances, create a single array with all combinations
if proc_type == 0:
    for exp_cont,exp in enumerate(IDs):
        cont_var_range_aux = []

        for ggex in g_ex_range[exp_cont]:
            for ggin in g_in_range[exp_cont]:
                cont_var_range_aux.append([ggex,ggin])

        cont_var_range.append(cont_var_range_aux)

# Cluster of results
cluster_0_rate_1 = [np.array([]),np.array([])]
cluster_1_rate_1 = [np.array([]),np.array([])]
cluster_2_rate_1 = [np.array([]),np.array([])]
cluster_0_rate_2 = [np.array([]),np.array([])]
cluster_1_rate_2 = [np.array([]),np.array([])]
cluster_2_rate_2 = [np.array([]),np.array([])]

# Average values
g_cluster_0 = []
g_cluster_1 = []
g_cluster_2 = []

print("Merging LFP traces...")

# Loop
for exp_cont,exp in enumerate(IDs):
    for v0 in ext_rates:

        for var in cont_var_range[exp_cont]:
        # for ggex in g_ex_range[exp_cont]:
        #     for ggin in g_in_range[exp_cont]:
            # Select only g values within canonical axes
            if proc_type == 0:
                if var[0]== 1. or var[1]== 1.:
                    # Print simulation case
                    var_ratio = -var[1]*Network_params["inh_exc_recurrent"]/(var[0]*Network_params["exc_exc_recurrent"])
                    # print("\nExt. rate = %s sp/s, g_ex = %s, g_in = %s, g = %s \n" % (v0,var[0],var[1],var_ratio))
                    filename = 'trial_'+str(0)+'_rate_'+str(v0)+'_gex_'+\
                                str(float("{0:.15f}".format(var[0])))+'_gin_'+str(float("{0:.15f}".format(var[1])))
                    continue_processing = True
                else:
                    continue_processing = False
            else:
                # Print simulation case
                var_ratio = var
                print("\nExt. rate = %s sp/s, variable = %s \n" % (v0,var))
                filename = 'trial_'+str(0)+'_rate_'+str(v0)+'_'+cont_var+'_'+\
                            str(float("{0:.15f}".format(var)))
                continue_processing = True

            # Load results
            if continue_processing:
                try:
                    with open('../results/'+exp+'/'+filename+".AMPA", "rb") as f:
                        AMPA_current = pickle.load(open('../results/'+exp+'/'+filename+".AMPA", "rb"),encoding='latin1')
                        GABA_current = pickle.load(open('../results/'+exp+'/'+filename+".GABA", "rb"),encoding='latin1')
                        t_sim = pickle.load(open('../results/'+exp+'/'+filename+".times", "rb"),encoding='latin1')
                        simstep = pickle.load(open('../results/'+exp+'/'+filename+".dt", "rb"),encoding='latin1')

                        AMPA_current = np.concatenate((AMPA_current,np.zeros(10)))
                        GABA_current = np.concatenate((GABA_current,np.zeros(10)))

                        if v0 == 1.5:
                            if var_ratio < 7.5 and len(cluster_0_rate_1[0])/100000 < 20:
                                g_cluster_0.append(var_ratio)
                                cluster_0_rate_1[0] = np.concatenate((cluster_0_rate_1[0],AMPA_current))
                                cluster_0_rate_1[1] = np.concatenate((cluster_0_rate_1[1],GABA_current))
                            elif var_ratio > 11 and len(cluster_2_rate_1[0])/100000 < 20:
                                g_cluster_2.append(var_ratio)
                                cluster_2_rate_1[0] = np.concatenate((cluster_2_rate_1[0],AMPA_current))
                                cluster_2_rate_1[1] = np.concatenate((cluster_2_rate_1[1],GABA_current))
                            else:
                                if len(cluster_1_rate_1[0])/100000 < 20:
                                    g_cluster_1.append(var_ratio)
                                    cluster_1_rate_1[0] = np.concatenate((cluster_1_rate_1[0],AMPA_current))
                                    cluster_1_rate_1[1] = np.concatenate((cluster_1_rate_1[1],GABA_current))
                        else:
                            if var_ratio < 7.5 and len(cluster_0_rate_2[0])/100000 < 20:
                                cluster_0_rate_2[0] = np.concatenate((cluster_0_rate_2[0],AMPA_current))
                                cluster_0_rate_2[1] = np.concatenate((cluster_0_rate_2[1],GABA_current))
                            elif var_ratio > 11 and len(cluster_2_rate_2[0])/100000 < 20:
                                cluster_2_rate_2[0] = np.concatenate((cluster_2_rate_2[0],AMPA_current))
                                cluster_2_rate_2[1] = np.concatenate((cluster_2_rate_2[1],GABA_current))
                            else:
                                if len(cluster_1_rate_2[0])/100000 < 20:
                                    cluster_1_rate_2[0] = np.concatenate((cluster_1_rate_2[0],AMPA_current))
                                    cluster_1_rate_2[1] = np.concatenate((cluster_1_rate_2[1],GABA_current))

                except Exception as e:
                    # print(e)
                    None

# Print number of samples
print("Rate 1.5 sp/s: %s, %s, %s" % (len(cluster_0_rate_1[0])/100000,len(cluster_1_rate_1[0])/100000,len(cluster_2_rate_1[0])/100000))
print("Rate 2 sp/s: %s, %s, %s" % (len(cluster_0_rate_2[0])/100000,len(cluster_1_rate_2[0])/100000,len(cluster_2_rate_2[0])/100000))
print("mean(g) = %s,%s,%s" % (np.mean(g_cluster_0),np.mean(g_cluster_1),np.mean(g_cluster_2)))

#######################################################################
# Save data
# os.mkdir('../results/merged_files')

# A) 1.5 sp/s
filename = 'trial_'+str(0)+'_rate_'+str(1.5)+'_gex_'+str(1.78289473684)+'_gin_'+str(1.0)
# AMPA and GABA
tools.saveData("merged_files",filename,".AMPA",cluster_0_rate_1[0])
tools.saveData("merged_files",filename,".GABA",cluster_0_rate_1[1])
# Save time array
tools.saveData("merged_files",filename,".times",t_sim)
# Save time step
tools.saveData("merged_files",filename,".dt",simstep)

filename = 'trial_'+str(0)+'_rate_'+str(1.5)+'_gex_'+str(1.20394736842)+'_gin_'+str(1.0)
# AMPA and GABA
tools.saveData("merged_files",filename,".AMPA",cluster_1_rate_1[0])
tools.saveData("merged_files",filename,".GABA",cluster_1_rate_1[1])
# Save time array
tools.saveData("merged_files",filename,".times",t_sim)
# Save time step
tools.saveData("merged_files",filename,".dt",simstep)

filename = 'trial_'+str(0)+'_rate_'+str(1.5)+'_gex_'+str(0.842105263158)+'_gin_'+str(1.0)
# AMPA and GABA
tools.saveData("merged_files",filename,".AMPA",cluster_2_rate_1[0])
tools.saveData("merged_files",filename,".GABA",cluster_2_rate_1[1])
# Save time array
tools.saveData("merged_files",filename,".times",t_sim)
# Save time step
tools.saveData("merged_files",filename,".dt",simstep)
########################################################################
# B) 2 sp/s
filename = 'trial_'+str(0)+'_rate_'+str(2.0)+'_gex_'+str(1.78289473684)+'_gin_'+str(1.0)
# AMPA and GABA
tools.saveData("merged_files",filename,".AMPA",cluster_0_rate_2[0])
tools.saveData("merged_files",filename,".GABA",cluster_0_rate_2[1])
# Save time array
tools.saveData("merged_files",filename,".times",t_sim)
# Save time step
tools.saveData("merged_files",filename,".dt",simstep)

filename = 'trial_'+str(0)+'_rate_'+str(2.0)+'_gex_'+str(1.20394736842)+'_gin_'+str(1.0)
# AMPA and GABA
tools.saveData("merged_files",filename,".AMPA",cluster_1_rate_2[0])
tools.saveData("merged_files",filename,".GABA",cluster_1_rate_2[1])
# Save time array
tools.saveData("merged_files",filename,".times",t_sim)
# Save time step
tools.saveData("merged_files",filename,".dt",simstep)

filename = 'trial_'+str(0)+'_rate_'+str(2.0)+'_gex_'+str(0.842105263158)+'_gin_'+str(1.0)
# AMPA and GABA
tools.saveData("merged_files",filename,".AMPA",cluster_2_rate_2[0])
tools.saveData("merged_files",filename,".GABA",cluster_2_rate_2[1])
# Save time array
tools.saveData("merged_files",filename,".times",t_sim)
# Save time step
tools.saveData("merged_files",filename,".dt",simstep)
