import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
#!#from mdtraj.reporters import HDF5Reporter 
import numpy as np
import sys

###############################################################################

def applyHarmonicPositionalRestraints(system, forceConstantInKcalPerMolePerAngSquared,
                                      inpcrd, indexOfAtomsToBeModified):
    """ This is essentially mimicking AMBER's restraint_wt"""

    forceConstant = u.Quantity(value=forceConstantInKcalPerMolePerAngSquared,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))

    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")

    force.addGlobalParameter("k",
       forceConstant.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))

    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for i in indexOfAtomsToBeModified:
        force.addParticle(i, inpcrd.getPositions()[i])

    system.addForce(force)

   
###############################################################################    
#################### Sym parameters ###########################################
###############################################################################

#!#prmtopName = 'PRMTOP'
#!#inpcrdName = 'INPCRD'

temperature = float(sys.argv[1])*u.kelvin
checkpoint_in_file = sys.argv[2]
csv_out_file = sys.argv[3]
dat_out_file = sys.argv[4]
pdb_out_file = sys.argv[5]
#!#traj_out_file = sys.argv[6]

#!#keyInteraction_ind_mol2 = [XXX, YYY]
keyInteraction = [keyInteraction_ind_mol2[0]-1, keyInteraction_ind_mol2[1]-1]

#!#spring_k = XX * u.kilocalorie/(u.mole * u.angstrom * u.angstrom)
#!#dist_in = XXX * u.angstrom # in angstrom 
#!#dist_fin = XXX * u.angstrom # in angstrom

steps_per_move = 200

velocity = 0.00001 * u.angstrom

# Platform definition

platform = mm.Platform_getPlatformByName("CUDA")
platformProperties = {}
platformProperties['CudaPrecision'] = 'mixed'
platformProperties['CudaDeviceIndex'] = '0'


# Read files

prmtop = app.AmberPrmtopFile(prmtopName)
inpcrd = app.AmberInpcrdFile(inpcrdName)


# Get indexes of heavy atoms in chunk

Chunk_Heavy_Atoms = []
for atom_i in prmtop.topology.atoms():
    if atom_i.residue.name not in ('HOH', 'WAT', 'IP3', 'LIG', 'Cs+', 'K+', 'Rb+', 'Li+', 'Na+', 'IP', 'Cl-', 'IM', 'IB') and atom_i.name[0] != 'H':
        Chunk_Heavy_Atoms.append(atom_i.index)
        

# Setting System

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds, hydrogenMass=None)


# Apply force on all havy atoms of chunk

applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)


# Integrator

integrator = mm.LangevinIntegrator(temperature, 4/u.picosecond, 0.002*u.picosecond)


# Setting Simulation object and loading the checkpoint

simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.loadCheckpoint(checkpoint_in_file)


# Get positions

positions = simulation.context.getState(getPositions=True).getPositions()


# SMD force definition

pullforce = mm.CustomExternalForce('k_sp*0.5*(R-R0)^2; \
                                   R = sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);')

pullforce.addPerParticleParameter('k_sp')

pullforce.addGlobalParameter("x0", 0.0 * u.nanometer)
pullforce.addGlobalParameter("y0", 0.0 * u.nanometer)
pullforce.addGlobalParameter("z0", 0.0 * u.nanometer)
pullforce.addGlobalParameter("R0", 0.0 * u.nanometer)

pullforce.addParticle(keyInteraction[1], [spring_k])

system.addForce(pullforce)


# Redefine integrator and simulation, and load checkpoint with new-updated system

integrator = mm.LangevinIntegrator(temperature, 4/u.picosecond, 0.002*u.picosecond)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.loadCheckpoint(checkpoint_in_file)


# Initializing energy

work_val_old = u.Quantity(value=0, unit=u.kilocalorie/u.mole)


# Number of big steps and pull distance 

steps = int(round((dist_fin  - dist_in) / velocity) / steps_per_move)
pull_distance = velocity * steps_per_move


# Reporters and duck.dat file

simulation.reporters.append(app.StateDataReporter(csv_out_file, steps_per_move, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, density=True, progress=True, totalSteps = steps_per_move * steps, speed=True))
#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))
f=open(dat_out_file,'w')


# Production in N steps with the update every 200 steps (2 pm) 

for i in range(steps):    

    # Get current state tu update the system

    state = simulation.context.getState(getPositions=True)
    pos_keyInt = state.getPositions()
    keyInteraction_pos = [pos_keyInt[keyInteraction[0]], pos_keyInt[keyInteraction[1]]]

    # Get radius of starting point and end point

    R_val = (dist_in + float(i+1) * pull_distance)
    R_val_start = (dist_in + float(i) * pull_distance)

    # Get distance of main interaction

    keyInteraction_dist = np.linalg.norm(keyInteraction_pos[0]-keyInteraction_pos[1])

    # Updated system

    simulation.context.setParameter('x0', keyInteraction_pos[0][0])
    simulation.context.setParameter('y0', keyInteraction_pos[0][1])
    simulation.context.setParameter('z0', keyInteraction_pos[0][2])
    simulation.context.setParameter('R0', R_val)

    # Calculate force F = -k * x

    force_val = -spring_k * (keyInteraction_dist - R_val)

    # Make step

    simulation.step(steps_per_move)

    # Calculate work for difference in potential energy in tranzition
    # W = EK_end - EK_start 
    # EK = 0.5 * k * x^2

    spr_energy_end = 0.5 * spring_k * (keyInteraction_dist - R_val)**2
    spr_energy_start = 0.5 * spring_k * (keyInteraction_dist - R_val_start)**2
    
    work_val = work_val_old + spr_energy_end - spr_energy_start

    work_val_old = work_val

    # Write duck.dat file 

    f.write(str(i)+' '+str(R_val)+' '+str(keyInteraction_dist)+' '+str(force_val)+' '+str(work_val)+'\n')
    
    
f.close()

# Save state in PDB file

positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out_file, 'w'))
