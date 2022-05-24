# Please note that this programm computes velocities v and not momenta p. However, in the publication we used
# a model system with mass m=1 and consequently p=v.
# Also: This script serves as an example to generate the Boltzmann and Maxwell-Boltzmann distribution for ABOBA.
# To compute the distributions for BAOAB, BAOA and GSD, change the argument in sim.compute_Langevin_dynamics()


import numpy as np
import matplotlib.pyplot as plt
from inspect import getsource
import sys

#############################################################################################
###   M D _ C L A S S   M O D E L   S Y S T E M   ###########################################
#############################################################################################

class Simulation_1D_model_system:

    def __init__(self, total_steps, potential, derivative_potential, time_step, collision_rate, mass, temperature, kboltz=0.008314, initial_position=0, initial_velocity=0):
        self.potential = potential
        self.derivative_potential = derivative_potential
        self.time_step = time_step
        self.collision_rate = collision_rate
        self.mass = mass
        self.temperature = temperature
        self.total_steps = total_steps
        self.kboltz = kboltz
        self.initial_position = initial_position
        self.initial_velocity = initial_velocity
        self.total_time = self.time_step * self.total_steps
        self._c1 = np.exp(-self.collision_rate*self.time_step)
        self._c2 = self.time_step/2
        self._c3 = self._c2/self.mass
        self._c4 = (self.kboltz*self.temperature)/self.mass
        self._c5 = np.exp(-(self.collision_rate*self.time_step)/2)

    def write_log_file(self, filepath):
        f= open(filepath + "traj_LOG.txt","w+")
        f.write("\n******* USED PARAMETERS ******** \n\n")
        f.write("time step = %s \n" %self.time_step)
        f.write("collision rate = %s \n" %self.collision_rate)
        f.write("mass = %s \n" %self.mass)
        f.write("temperature = %s \n" %self.temperature)
        f.write("total number of steps = %s \n" %self.total_steps)
        f.write("Boltzmann constant = %s \n" %self.kboltz)
        f.write("initial position = %s \n" %self.initial_position)
        f.write("initial velocity = %s \n" %self.initial_velocity)
        f.write("%s" %getsource(self.potential))
        f.write("%s" %getsource(self.derivative_potential))
        f.close()

    def draw_random_number_sequence(self, mean=0, variance=1, length=None):
        if length == None:
            eta = np.random.normal(mean, variance, (self.total_steps))
        else:
            eta = np.random.normal(mean, variance, (length))
        return eta

    def compute_Langevin_dynamics(self, random_numbers=np.array([]), integrator=''):
        position = np.zeros((self.total_steps + 1))
        position[0] = self.initial_position
        velocity = np.zeros((self.total_steps + 1))
        velocity[0] = self.initial_velocity
        #---- check for random number sequence, else generate it  --------------------------
        if random_numbers.size == 0:
            eta = self.draw_random_number_sequence()
        else:
            eta = random_numbers
        #--- raise error if integrator not given at all ------------------------------------
        if integrator == '':
            sys.exit('Error in compute_Langevin_dynamics: Please specifiy integration scheme!')
        #--- BAOAB integrator -------------------------------------------------------------
        elif integrator == 'BAOAB':
            for k in range(0, self.total_steps):
                velocity[k+1] = velocity[k] - self._c3*self.derivative_potential(position[k])
                position[k+1] = position[k] + self._c2*velocity[k+1]
                velocity[k+1] = self._c1*velocity[k+1] + np.sqrt(self._c4*(1 - self._c1**2))*eta[k]
                position[k+1] = position[k+1] + self._c2*velocity[k+1]
                velocity[k+1] = velocity[k+1] - self._c3*self.derivative_potential(position[k+1])
        #--- ABOBA integrator -------------------------------------------------------------
        elif integrator == 'ABOBA':
            for k in range(0, self.total_steps):
                position[k+1] = position[k] + self._c2*velocity[k]
                velocity[k+1] = velocity[k] - self._c3*self.derivative_potential(position[k+1])
                velocity[k+1] = self._c1*velocity[k+1] + np.sqrt(self._c4*(1 - self._c1*self._c1))*eta[k]
                velocity[k+1] = velocity[k+1] - self._c3*self.derivative_potential(position[k+1])
                position[k+1] = position[k+1] + self._c2*velocity[k+1]
        #--- GSD ---------------------------------------------------------------------------
        elif integrator == 'GSD':
            for k in range(0, self.total_steps):
                velocity[k+1] = velocity[k] - (self.time_step/self.mass)*self.derivative_potential(position[k])
                dv = -(1 - self._c1)*velocity[k+1] + np.sqrt(self._c4*(1 - self._c1*self._c1))*eta[k]
                position[k+1] = position[k] + (velocity[k+1] + 0.5*dv)*self.time_step
                velocity[k+1] = velocity[k+1] + dv
        #--- BAOA --------------------------------------------------------------------------
        elif integrator == 'BAOA':
            for k in range(0, self.total_steps):    
                velocity[k+1] = velocity[k] - (self.time_step/self.mass)*self.derivative_potential(position[k])
                position[k+1] = position[k] + self._c2*velocity[k+1]
                velocity[k+1] = self._c1*velocity[k+1] + np.sqrt(self._c4*(1 - self._c1*self._c1))*eta[k]
                position[k+1] = position[k+1] + self._c2*velocity[k+1]
        #--- raise error if integrator not recognized or implemented -----------------------
        else:
            sys.exit('Error in compute_Langevin_dynamics: Integrator not implemented!')
        return position, velocity



#############################################################################################
###   P A R A M E T E R S  +  I N P U T   ###################################################
#############################################################################################

#---- Simulation parameters -------------------------------------------------------
nsteps = 5*10**7  # total number of steps
V = lambda x: (x**2 - 1)**2 + x  # potential
der_V = lambda x: 4*x*(x**2 - 1) + 1 # derivative of potential
dt = 0.25   # time step
xi = 1      # collision rate
m = 1       # mass
T = 1       # temperature
kboltz = 1  # Boltzmann constant

#--- Distribution parameters ------------------------------------------------------
Boltzmann_nbins = 100   # number of bins
Boltzmann_min = -2      # minimal value position
Boltzmann_max = 2       # maximal value position
Maxwell_nbins = 100     # number of bins
Maxwell_min = -5        # minimal value velocity
Maxwell_max = 5         # maximal value velocity   



#############################################################################################
###   G E N E R A T E   T R A J E C T O R Y   ###############################################
#############################################################################################

# initialize simulation object
sim = Simulation_1D_model_system(total_steps = nsteps,
         potential=V, 
         derivative_potential=der_V,
         time_step=dt,
         collision_rate=xi,
         mass=m,
         temperature=T,
         kboltz=kboltz)

# write simulation LOG-file
sim.write_log_file("./") 
# write distribution LOG-file
f= open("distribution_LOG.txt","w+")
f.write("\n******* USED PARAMETERS ******** \n\n")
f.write("Boltzmann distribution \n")
f.write("------------------------------------------------------- \n")
f.write("minimal position = %s \n" %Boltzmann_min)
f.write("maximal position = %s \n" %Boltzmann_max)
f.write("number of bins = %s \n" %Boltzmann_nbins)
f.write("\nMaxwell-Boltzmann distribution \n")
f.write("------------------------------------------------------- \n")
f.write("minimal velocity = %s \n" %Maxwell_min)
f.write("maximal velocity = %s \n" %Maxwell_max)
f.write("number of bins = %s \n" %Maxwell_nbins)
f.close()

# generate trajectory for ABOBA
# For BAOAB, BAOA or GSD, change the integrator string.
x, v = sim.compute_Langevin_dynamics(integrator="ABOBA")



#############################################################################################
###   A N A L Y S I S   #####################################################################
#############################################################################################

# ---- Compute distributions --------------------------------------------------

# prepare arrays to save results
results_Boltzmann = np.zeros((2, Boltzmann_nbins))
results_Maxwell = np.zeros((2, Maxwell_nbins))
    
# compute Boltzmann distribution (positions) 
results_Boltzmann[1], bins_x = np.histogram(x, bins=Boltzmann_nbins, range=(Boltzmann_min,Boltzmann_max), density=True)
results_Boltzmann[0] = bins_x[:-1] + (bins_x[-1] - bins_x[-2])/2

# compute Maxwell-Boltzmann distribution (velocities)
results_Maxwell[1], bins_v = np.histogram(v, bins=Maxwell_nbins, range=(Maxwell_min,Maxwell_max), density=True)
results_Maxwell[0] = bins_v[:-1] + (bins_v[-1] - bins_v[-2])/2



#############################################################################################
###   S A V E S   ###########################################################################
#############################################################################################

# change name for BAOAB, BAOA or GSD
#np.save("results/distribution_positions_ABOBA.npy", results_Boltzmann)
#np.save("results/distribution_velocities_ABOBA.npy", results_Maxwell)













