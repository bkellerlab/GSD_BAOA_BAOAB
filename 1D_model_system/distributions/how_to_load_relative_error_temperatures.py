import numpy as np
import matplotlib.pyplot as plt



#############################################################################################
###   L O A D   #############################################################################
#############################################################################################

integrators = ["ABOBA", "BAOAB", "GSD", "BAOA"]
Tkin = np.load('results/relative_error_kinetic_temp.npy', allow_pickle=True).item()
Tpos = np.load('results/relative_error_configurational_temp.npy', allow_pickle=True).item()


#############################################################################################
###   IP L O T   ############################################################################
#############################################################################################

# example relative error in kinetic temperature

for ii in range(0, len(integrators)):
    plt.plot(Tkin[integrators[ii]][0], Tkin[integrators[ii]][1], label=integrators[ii])
plt.yscale('log')
plt.legend()






































