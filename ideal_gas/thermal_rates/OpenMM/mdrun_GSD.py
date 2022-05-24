# This script is for the GSD integrator and a collision rate of 1 1/ps
# change the value of the variable xi for a different collision rate

import openmmtools as tools
import simtk.openmm as openmm
from simtk import unit
from time import gmtime, strftime
import numpy as np


##########################################################################################################
###  I N P U T   A N D   S E T - U P  ####################################################################
##########################################################################################################

# ideal gas in box parameters
nparticles = 3200
mass = 72*unit.amu # mass of each particle
box_volume = 87.19095183464606*(unit.nanometers**3)

# equilibration parameters
T_equi = 350*unit.kelvin  # equilibration temperature
equisteps = 1000  # total number of equilibration steps

# simulation parameters
T_sim = 320*unit.kelvin  # production run temperature
xi = 1*(unit.picoseconds)**(-1)  # collision rate
dt = 2*unit.femtoseconds  # time step
nsteps = 25000  # total number of steps for production run
nstxout = 1  # write-out frequency in steps
kboltz = 0.008314*unit.kilojoule/(unit.mole*unit.kelvin) # Boltzmann constant

# create system
ideal_gas = tools.testsystems.IdealGas(nparticles=nparticles, mass=mass, temperature=T_equi, volume=box_volume)
system, initial_positions, topology = ideal_gas.system, ideal_gas.positions, ideal_gas.topology

# get hardware stuff done
platform = openmm.Platform.getPlatformByName('CPU')

# write LOG-file
log = open("LOG.txt", 'w')
log.write('\nSystem parameters: \n')
log.write('--------------------------\n')
log.write('number of particles: ' + str(nparticles) + "\n")
log.write("mass per particle: " + str(mass) + " (a.u.) \n" )
log.write("box volume: " + str(box_volume) + "\n" )
log.write('\n\nEquilibration parameters: \n')
log.write('--------------------------\n')
log.write('temperature: ' + str(T_equi) + "\n")
log.write('steps: ' + str(equisteps) + " steps \n")
log.write('time: ' + str(equisteps*dt) + "\n")
log.write('\n\nSimulation parameters: \n')
log.write('--------------------------\n')
log.write('temperature: ' + str(T_sim) + "\n")
log.write('collision rate: ' + str(xi) + "\n")
log.write('Boltzmann constant: ' + str(kboltz) + "\n")
log.write('time step: ' + str(dt) + "\n")
log.write('steps: ' + str(nsteps) + " steps \n")
log.write('time: ' + str(nsteps*dt) + "\n")
log.write('write-out frequency: ' + str(nstxout) + " steps \n")
log.write("Simulation start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n" )
log.close();



##########################################################################################################
###  G S D   I N T E G R A T O R   #######################################################################
##########################################################################################################

# note that there are no constraints because ideal gas --> no constraints handling necessary

# integrator constants
c1 = np.exp(-xi*dt)
c2 = np.sqrt(kboltz*T_equi*(1 - np.exp(-2*xi*dt)))

# initialize integrator 
integrator = openmm.openmm.CustomIntegrator(dt)

# define per degree of freedom (PerDof) and global variables
integrator.addGlobalVariable("c1", c1)
integrator.addGlobalVariable("c2", c2)
integrator.addGlobalVariable("dt", dt)
integrator.addPerDofVariable("eta", 0)
integrator.addPerDofVariable("dv", 0)

#------ integrator equations --------------------------------------------------

# draw random number for upcoming step
integrator.addComputePerDof("eta","gaussian")

# GSD step 1: intermediate velocities
integrator.addComputePerDof("v", "v + f*(dt/m)")
# GSD step 2: velocity increment
integrator.addComputePerDof("dv", "-v + c1*v + c2*sqrt(1/m)*eta")
# GSD step 3: position update
integrator.addComputePerDof("x", "x + (v + 0.5*dv)*dt")
# GSD step 4: velocitiy update
integrator.addComputePerDof("v", "v + dv")



##########################################################################################################
###  E Q U I L I B R A T I O N   #########################################################################
##########################################################################################################

# set-up simulation
simulation = openmm.app.simulation.Simulation(topology, system, integrator, platform)
simulation.context.setPositions(initial_positions)

# minimization
print('\n\n*** Minimizing ...')
simulation.minimizeEnergy()

# equilibration
print('\n\n*** Equilibrating...')
simulation.context.setVelocitiesToTemperature(T_equi)
simulation.step(equisteps)

# get positions and velocities from last equilibration frame
state = simulation.context.getState(getPositions=True, getVelocities=True)
equi_positions = state.getPositions()
equi_velocities = state.getVelocities()



##########################################################################################################
###  S I M U L A T I O N   A N D   O U T P U T   #########################################################
##########################################################################################################

print('\n\n*** Simulating ...')
# change temperature in integrator for simulation
c2 = np.sqrt(kboltz*T_sim*(1 - np.exp(-2*xi*dt)))
integrator.setGlobalVariableByName("c2", c2)

# use positions and velocities of last equilibration frame as starting point for production run
simulation.context.setPositions(equi_positions)
simulation.context.setVelocities(equi_velocities)

# save trajectory
#simulation.reporters.append(openmm.app.pdbreporter.PDBReporter("trajectory.pdb", nstxout))
# save instantaneous temperature and others
simulation.reporters.append(openmm.app.statedatareporter.StateDataReporter("data.txt", nstxout, step=True, kineticEnergy=True, temperature=True, separator=' '))

# repeat procedure for nsteps
simulation.step(nsteps)

# add total calculation time to LOG-file
log = open("LOG.txt", 'a')
log.write("Simulation end: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) )
log.close();

# end
print('\n****** Simulation Complete *****************************\n\n')


