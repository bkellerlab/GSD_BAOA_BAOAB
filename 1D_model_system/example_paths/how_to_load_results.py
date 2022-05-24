import numpy as np
import matplotlib.pyplot as plt



#############################################################################################
###   L O A D   D A T A   ###################################################################
#############################################################################################

integrators =["ABOBA", "BAOAB", "GSD", "BAOA"]

x = np.load('results/positions.npy', allow_pickle=True).item()
v = np.load('results/velocities.npy', allow_pickle=True).item()
dev_x = np.load('results/deviation_positions.npy', allow_pickle=True).item()
dev_v = np.load('results/deviation_velocities.npy', allow_pickle=True).item()
pathlength = np.load('results/pathlength.npy')


#############################################################################################
###   P L O T   D A T A   ###################################################################
#############################################################################################

# example: plot positions -----------------------------------------------------

fig1 = plt.figure()
for ii in range(0, len(integrators)):
    plt.plot(pathlength, x[integrators[ii]], label=integrators[ii])
plt.legend()

# equivalent for vleocities and deviations in positions and velocities.


























