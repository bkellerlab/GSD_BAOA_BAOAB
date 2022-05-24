# This script computes the average kinetic and average configurational temperature for ABOBA, BAOAB,
# BAOA and GSD and a fixed time step (here dt=0.2). You can execute n independent experiments per integrator 
# with this script. 

import numpy as np
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

    def compute_average_temperature(self, nexperiment, integrator=''):
        # compute Tkin and Tconf for nexperiment experiments in parallel
        #--- initialize starting conditions ------------------------------------------------
        position = np.ones((nexperiment))*self.initial_position
        velocity = np.ones((nexperiment))*self.initial_velocity
        Tkin = self.mass*velocity*velocity
        Tpos = position*self.derivative_potential(position)
        #--- raise error if integrator not given at all ------------------------------------
        if integrator == '':
            sys.exit('Error in compute_Langevin_dynamics: Please specifiy integration scheme!')
        #--- ABOBA integrator -------------------------------------------------------------
        elif integrator == 'ABOBA':
            for k in range(0, self.total_steps):
                eta = np.random.normal(loc=0, scale=1.0, size=nexperiment)
                position = position + self._c2*velocity
                velocity = velocity - self._c3*self.derivative_potential(position)
                velocity = self._c1*velocity + np.sqrt(self._c4*(1 - self._c1*self._c1))*eta
                velocity = velocity - self._c3*self.derivative_potential(position)
                position = position + self._c2*velocity
                # compute average kinetic temperature
                Tkin = Tkin + self.mass*velocity*velocity
                # compute average configurational temperature
                Tpos = Tpos + position*self.derivative_potential(position)
        #--- BAOAB integrator -------------------------------------------------------------
        elif integrator == 'BAOAB':
            for k in range(0, self.total_steps):
                eta = np.random.normal(loc=0, scale=1.0, size=nexperiment)
                velocity = velocity - self._c3*self.derivative_potential(position)
                position = position + self._c2*velocity
                velocity = self._c1*velocity + np.sqrt(self._c4*(1 - self._c1**2))*eta
                position = position + self._c2*velocity
                velocity = velocity - self._c3*self.derivative_potential(position)
                # compute average kinetic temperature
                Tkin = Tkin + self.mass*velocity*velocity
                # compute average configurational temperature
                Tpos = Tpos + position*self.derivative_potential(position)                             
        #--- GSD ---------------------------------------------------------------------------
        elif integrator == 'GSD':
            for k in range(0, self.total_steps):
                eta = np.random.normal(loc=0, scale=1.0, size=nexperiment)
                velocity = velocity - (self.time_step/self.mass)*self.derivative_potential(position)
                dv = -(1 - self._c1)*velocity + np.sqrt(self._c4*(1 - self._c1*self._c1))*eta
                position = position + (velocity + 0.5*dv)*self.time_step
                velocity = velocity + dv
                # compute average kinetic temperature
                Tkin = Tkin + self.mass*velocity*velocity
                # compute average configurational temperature
                Tpos = Tpos + position*self.derivative_potential(position)
        #--- BAOA --------------------------------------------------------------------------
        elif integrator == 'BAOA':
            for k in range(0, self.total_steps): 
                eta = np.random.normal(loc=0, scale=1.0, size=nexperiment)
                velocity = velocity - (self.time_step/self.mass)*self.derivative_potential(position)
                position = position + self._c2*velocity
                velocity = self._c1*velocity + np.sqrt(self._c4*(1 - self._c1*self._c1))*eta
                position = position + self._c2*velocity
                # compute average kinetic temperature
                Tkin = Tkin + self.mass*velocity*velocity
                # compute average configurational temperature
                Tpos = Tpos + position*self.derivative_potential(position)                
        #--- raise error if integrator not recognized --------------------------------------
        else:
            sys.exit('Error in compute_Langevin_dynamics: Please specifiy integration scheme!')
        return Tkin/(self.total_steps+1), Tpos/(self.total_steps+1)  
    
    
#############################################################################################
###   P A R A M E T E R S  +  I N P U T   ###################################################
#############################################################################################

#---- Simulation parameters -------------------------------------------------------
nsteps = 10**6  # total number of steps
V = lambda x: (x**2 - 1)**2 + x  # potential
der_V = lambda x: 4*x*(x**2 - 1) + 1 # derivative of potential
dt = 0.2    # time step
xi = 1      # collision rate
m = 1       # mass
T = 1       # temperature
kboltz = 1  # Boltzmann constant

# integrators
integrators = ["ABOBA", "BAOAB", "GSD", "BAOA"]
# number of independent experiments per integrator
nexperiment = 500



#############################################################################################
###   C O M P U T A T I O N   ###############################################################
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

# write log file
sim.write_log_file("./") 
# add number of independent experiments to LOG-file
f= open("traj_LOG.txt","a")
f.write("number of experiments = %s \n" %nexperiment)
f.close()

# --- compute average kinetic and average configurational temperature from simulation -----

Tkin = {}
Tpos = {}
for integrator in integrators:
    Tkin[integrator], Tpos[integrator] = sim.compute_average_temperature(nexperiment, integrator=integrator)



#############################################################################################
###   S A V E S   ###########################################################################
#############################################################################################

# save temperatures
#np.save("results/average_configurational_temp_dt_020.npy", Tpos)
#np.save("results/average_kinetic_temp_dt_020.npy", Tkin)


