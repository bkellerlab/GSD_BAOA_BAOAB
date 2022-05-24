# Please add the path that points to the forcefield you want to use (we used forcefield from Gromacs library)
# to the variable forcefield_directory
#---------------------------------------------------------------------------------------------------------
# This script can be adapted for BAOAB and BAOA.
# For BAOAB: change splitting string in tools.integrators.LangevinIntegrator() to splitting='V R O R V'
# For BAOA: change splitting string in tools.integrators.LangevinIntegrator() to splitting='V R O R'

import openmmtools as tools
import simtk.openmm as openmm
from simtk import unit
from time import gmtime, strftime
import numpy as np

#%%

##########################################################################################################
###  I N P U T   A N D   S E T - U P  ####################################################################
##########################################################################################################

# paths
folder_input = "../input/"
forcefield_directory = 'insert path to forcefield you want to use'

# water box parameters (just for LOG-file, not used in simulation script because information is already inclduded in input files)
cutoff = 1*unit.nanometer
boxlength = 3.143*unit.nanometer
nmolecules = 1024

# equilibration parameters
equisteps = 2000

# simulation parameters
T = 300*unit.kelvin  # temperature
xi = 2*(unit.picoseconds)**(-1)  # collision rate
dt = 2*unit.femtoseconds  # time step
nsteps = 2500000  # total number of steps of production run
nstxout = 10      # write-out frequency in steps

#Input topology+coordinates files            
gro = openmm.app.gromacsgrofile.GromacsGroFile(folder_input + 'water_box.gro')
top = openmm.app.gromacstopfile.GromacsTopFile(folder_input + 'water_box.top', periodicBoxVectors=gro.getPeriodicBoxVectors(),includeDir=forcefield_directory)

#create system                         
system = top.createSystem(nonbondedMethod=openmm.app.PME, nonbondedCutoff=cutoff, constraints=openmm.app.HBonds)

# get hardware stuff done
platform = openmm.Platform.getPlatformByName('CPU')


##########################################################################################################
###  L O G - F I L E   ###################################################################################
##########################################################################################################

# write LOG-file
log = open("LOG.txt", 'w')
log.write('\nSystem parameters: \n')
log.write('--------------------------\n')
log.write("box length: " + str(boxlength) + "\n" )
log.write("non-bonded cutoff: " + str(cutoff) + "\n" )
log.write("number of water molecules: " + str(nmolecules) + "\n" )
log.write('\n\nEquilibration parameters: \n')
log.write('--------------------------\n')
log.write('steps: ' + str(equisteps) + " steps \n")
log.write('time: ' + str(equisteps*dt) + "\n")
log.write('\n\nSimulation parameters: \n')
log.write('--------------------------\n')
log.write('temperature: ' + str(T) + "\n")
log.write('collision rate: ' + str(xi) + "\n")
log.write('time step: ' + str(dt) + "\n")
log.write('steps: ' + str(nsteps) + " steps \n")
log.write('time: ' + str(nsteps*dt) + "\n")
log.write('write-out frequency: ' + str(nstxout) + " steps \n")
log.write("Simulation start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n" )
log.close();


##########################################################################################################
###  S E T - U P   S I M U L A T I O N   #################################################################
##########################################################################################################

# set-up integrator 
integrator = tools.integrators.LangevinIntegrator(temperature=T, collision_rate=xi, timestep=dt, splitting='R V O V R')
# set-up simulation
simulation = openmm.app.simulation.Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)



##########################################################################################################
###  E Q U I L I B R A T I O N   #########################################################################
##########################################################################################################

# minimization
print('\n\n*** Minimizing ...')
simulation.minimizeEnergy()

# equilibration
print('\n\n*** Equilibrating...')
simulation.context.setVelocitiesToTemperature(T)
simulation.step(equisteps)

##########################################################################################################
###  S I M U L A T I O N   A N D   O U T P U T   #########################################################
##########################################################################################################

print('\n\n*** Simulating ...')

# save trajectory (uncomment if trajectory should be written out)
#simulation.reporters.append(openmm.app.pdbreporter.PDBReporter("trajectory.pdb", nstxout))
# save instantaneous temperature and others
simulation.reporters.append(openmm.app.statedatareporter.StateDataReporter("data.txt", nstxout, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, separator=' '))

# repeat procedure for nsteps
simulation.step(nsteps)

# add total calculation time to LOG-file
log = open("LOG.txt", 'a')
log.write("Simulation end: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) )
log.close();

# end
print('\n****** Simulation Complete *****************************\n\n')







