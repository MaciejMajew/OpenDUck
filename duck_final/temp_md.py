import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
#!#from mdtraj.reporters import HDF5Reporter
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

#!#prmtopName = 'PRMTOP'
#!#inpcrdName = 'INPCRD'

checkpoint_in_file = sys.argv[1]
checkpoint_out_file = sys.argv[2]
csv_out_file = sys.argv[3]
pdb_out_file = sys.argv[4]
#!#traj_out_file = sys.argv[5]

#!#MD_len = X * u.nanosecond
sim_steps = round(MD_len / (0.002 * u.picosecond))

#!#keyInteraction_ind_mol2 = [XXX, YYY]
keyInteraction = [keyInteraction_ind_mol2[0]-1, keyInteraction_ind_mol2[1]-1]


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


# Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance

applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)
applyLigandChunkRestraint(system, 1.0, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)


# Integrator

integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)


# Setting Simulation object and loading the checkpoint

simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.loadCheckpoint(checkpoint_in_file)


# Simulation reporters

simulation.reporters.append(app.StateDataReporter(csv_out_file, 2000, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, density=True, progress=True, totalSteps=250000, speed=True))
#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))


# Production

simulation.step(sim_steps)


# Save state in checkpoint file and save coordinates in PDB file
            
state = simulation.context.getState(getPositions=True, getVelocities=True)
positions = state.getPositions()
velocities = state.getVelocities()

app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out_file, 'w'))

simulation.saveCheckpoint(checkpoint_out_file)
