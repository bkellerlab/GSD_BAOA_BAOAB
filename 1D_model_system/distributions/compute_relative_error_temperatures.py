# This script calculates the relative error for the kinetic and the configurational temperatues.
# Please generate the average kinetic and everage configurational temperature for n independent 
# experiments per integrator and different time steps dt first (use compute_temperatures.py), 
# before excuting this script.

import numpy as np
import matplotlib.pyplot as plt

#%%

#############################################################################################
###   I N P U T   ###########################################################################
#############################################################################################

integrators = ["ABOBA", "BAOAB", "GSD", "BAOA"]
# list with used labels for the saves that contain the average kinetic and configurationa temperatures
dts_labels = ["020", "022", "024", "026"]
# used time steps
dts = np.array([0.2, 0.22, 0.24, 0.26])
# reference temperature
Tref = 1



#############################################################################################
###   C O M P U T A T I O N   ###############################################################
#############################################################################################

# --- kinetic temperature ------------------------------
# compute the average from the n independent experiments per integrator
av_Tkin = {}
for integrator in integrators:
    # store time step and mean
    av_Tkin[integrator] = np.array([[], []])

for ii in range(0, len(dts)):
    Tkin = np.load("results/average_kinetic_temp_dt_" + dts_labels[ii] + ".npy", allow_pickle=True).item()
    for integrator in integrators:
        mean = np.mean(Tkin[integrator])
        # check if mean is a real number, if yes proceed
        if np.isnan(mean) == False:
            data = np.array([[dts[ii]], [abs(Tref - mean)]])
            av_Tkin[integrator] = np.append(av_Tkin[integrator], data, axis=1)
        

# --- configurational temperature ------------------------------
# compute the average from the n independent experiments per integrator
av_Tpos = {}
for integrator in integrators:
    # store time step and mean per integrator
    av_Tpos[integrator] = np.array([[], []])

for ii in range(0, len(dts)):
    Tpos = np.load("results/average_configurational_temp_dt_" + dts_labels[ii] + ".npy", allow_pickle=True).item()
    for integrator in integrators:
        mean = np.mean(Tpos[integrator])
        if np.isnan(mean) == False:
            data = np.array([[dts[ii]], [abs(Tref - mean)]])
            av_Tpos[integrator] = np.append(av_Tpos[integrator], data, axis=1)



#############################################################################################
###   S A V E S   ###########################################################################
#############################################################################################

#np.save('results/relative_error_kinetic_temp.npy', av_Tkin)
#np.save('results/relative_error_configurational_temp.npy', av_Tpos)   
    
    
    

