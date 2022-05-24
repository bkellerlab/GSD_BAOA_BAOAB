import numpy as np
import matplotlib.pyplot as plt



#############################################################################################
###   L O A D   #############################################################################
#############################################################################################

integrators = ["ABOBA", "BAOAB", "GSD", "BAOA"]

# example: potential energy at T=300 K
distributions = np.load('Epot_distribution_T_300.npy', allow_pickle=True).item()    



#############################################################################################
###   I N P U T   ###########################################################################
#############################################################################################

for ii in range(0, len(integrators)):
    plt.plot(distributions[integrators[ii]][0], distributions[integrators[ii]][1], label=integrators[ii])
plt.legend()






































