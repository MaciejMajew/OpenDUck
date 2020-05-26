import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u

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

def applyLigandChunkRestraint(system, k2, k3, R2, R3, R4, indexOfAffectedAtoms):

    forceConstant_k2 = u.Quantity(value=k2,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))
    forceConstant_k3 = u.Quantity(value=k3,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))

    restraint_force = mm.CustomBondForce('step(R2 - r) * f1 + step(r - R3) * select( step(r - R4), f3, f2);'
                                         'f1 = k2 * (r - R2)^2;'
                                         'f2 = k3 * (r - R3)^2;'
                                         'f3 = k3 * (R4 - R3) * (2 * r - R4 - R3)')

    restraint_force.addGlobalParameter("k2",
       forceConstant_k2.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))
    restraint_force.addGlobalParameter("k3",
       forceConstant_k3.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))

    restraint_force.addGlobalParameter("R2", R2)
    restraint_force.addGlobalParameter("R3", R3)
    restraint_force.addGlobalParameter("R4", R4)

    restraint_force.addBond(indexOfAffectedAtoms[0], indexOfAffectedAtoms[1])

    system.addForce(restraint_force)


###############################################################################
#################### Sym parameters ###########################################
###############################################################################

#!#keyInteraction_ind_mol2 = [XXX, YYY]
keyInteraction = [keyInteraction_ind_mol2[0]-1, keyInteraction_ind_mol2[1]-1]

# Platform definition

platform = mm.Platform_getPlatformByName("CUDA")
platformProperties = {}
platformProperties['CudaPrecision'] = 'mixed'
platformProperties['CudaDeviceIndex'] = '0'


# Read files

#!#prmtopName = 'PRMTOP'
#!#inpcrdName = 'INPCRD'

prmtop = app.AmberPrmtopFile(prmtopName)
inpcrd = app.AmberInpcrdFile(inpcrdName)


##################
##################
#  Minimisation  #
##################
##################

print('Minimising...')

# Define system

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom)

# Get indexes of heavy atoms in chunk

Chunk_Heavy_Atoms = []
for atom_i in prmtop.topology.atoms():
    if atom_i.residue.name not in ('HOH', 'WAT', 'IP3', 'LIG', 'Cs+', 'K+', 'Rb+', 'Li+', 'Na+', 'IP', 'Cl-', 'IM', 'IB') and atom_i.name[0] != 'H':
        Chunk_Heavy_Atoms.append(atom_i.index)

# Apply force on all havy atoms of chunk

applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)

# Integrator

integrator = mm.VerletIntegrator(1*u.femtosecond)

# Define Simulation

simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(inpcrd.positions)

# Minimizing energy

simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions

positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))


##########################
##########################
# Equlibration - heating #
##########################
##########################
#new minimised positions, however using old restraints

# Define new system

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds, hydrogenMass=None)

# Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance

applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)
applyLigandChunkRestraint(system, 1.0, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)

# Intergator 

integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)

# Define Simulation

simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions) #changing coordintes to minimized   

# Reporters

simulation.reporters.append(app.StateDataReporter("heating.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))

# Heating the system

simulation.step(50000) # 0.1 ns 

# Save the positions and velocities

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

app.PDBFile.writeFile(simulation.topology, positions, open('heating_final.pdb', 'w'))

#clear reporters

simulation.reporters = []

##########################
##########################
# Equlibration - density #
##########################
##########################

# Add barostat to the system

system.addForce(mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin))

# Integrator 

integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)

# Define simulation

simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

# Reporters

simulation.reporters.append(app.StateDataReporter("density.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))

# Correcting the density

simulation.step(50000) # 0.1 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
app.PDBFile.writeFile(simulation.topology, positions, open('density_final.pdb', 'w'))

#saving simulation stage
simulation.saveCheckpoint('equil.chk')
