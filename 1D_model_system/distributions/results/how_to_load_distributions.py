import numpy as np
import matplotlib.pyplot as plt



#############################################################################################
###   L O A D   #############################################################################
#############################################################################################

integrators = ["ABOBA", "BAOAB", "GSD", "BAOA"]

# example: positions distributions
distributions = {}
for integrator in integrators:
    distributions[integrator] = np.load('distribution_positions_' + integrator + '.npy')    



#############################################################################################
###   I N P U T   ###########################################################################
#############################################################################################

for ii in range(0, len(integrators)):
    plt.plot(distributions[integrators[ii]][0], distributions[integrators[ii]][1], label=integrators[ii])
plt.legend()






































