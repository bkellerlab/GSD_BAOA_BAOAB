import numpy as np
import matplotlib.pyplot as plt



##########################################################################################################
###  L O A D   D A T A   #################################################################################
##########################################################################################################

k = np.load('thermal_rates.npy', allow_pickle=True).item()
xis = np.load('collision_rates.npy')
integrators = ["ABOBA", "BAOAB", "GSD", "BAOA"]

# analytic thermal rate for ideal gas
dt = 0.002  # used time step
k_ana = lambda xi: (1 - np.exp(-2*xi*dt))/dt




##########################################################################################################
###  P L O T   ###########################################################################################
##########################################################################################################

for ii in range(0, len(integrators)):
    plt.plot(xis, k_ana(xis), ':', color='r', linewidth=1)  
    plt.errorbar(xis, k[integrators[ii]][0], yerr=k[integrators[ii]][1], fmt='o', capsize=3, label=integrators[ii])
plt.legend()













