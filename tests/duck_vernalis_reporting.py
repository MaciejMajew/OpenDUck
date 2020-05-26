#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import mdtraj as md
from mdtraj.reporters import HDF5Reporter
import numpy as np


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
    
def SMD(system, prmtop, platform, platformProperties, temperature, positions, velocities, keyInteraction, spring_k, dist_in, dist_fin, SMD_num, save_step, move_force_step):

    # See page 456 of http://ambermd.org/doc12/Amber17.pdf
    pullforce = mm.CustomExternalForce('k_sp*0.5*(dx^2+dy^2+dz^2); \
                                       dx=x-(x0+displace_x); \
                                       dy=y-(y0+displace_y); \
                                       dz=z-(z0+displace_z);')

    pullforce.addPerParticleParameter('k_sp')
    pullforce.addPerParticleParameter('x0')
    pullforce.addPerParticleParameter('y0')
    pullforce.addPerParticleParameter('z0')

    pullforce.addGlobalParameter("displace_x", 0.0 * u.nanometer)
    pullforce.addGlobalParameter("displace_y", 0.0 * u.nanometer)
    pullforce.addGlobalParameter("displace_z", 0.0 * u.nanometer)

    keyInteraction_pos = [positions[keyInteraction[0]], positions[keyInteraction[1]]]
    keyInteraction_dist = np.linalg.norm(keyInteraction_pos[0] - keyInteraction_pos[1])
    keyInteraction_vect = (keyInteraction_pos[1] - keyInteraction_pos[0]) / keyInteraction_dist
    keyInteraction_vect = u.Quantity(value=keyInteraction_vect, unit=u.nanometers)
    pullto = keyInteraction_pos[0] + 0.25 * keyInteraction_vect
    pullto_old = pullto

    pullforce.addParticle(keyInteraction[1], [spring_k, pullto[0], pullto[1], pullto[2] ])
    system.addForce(pullforce)

    integrator = mm.LangevinIntegrator(temperature, 4/u.picosecond, 0.002*u.picosecond)

    simulation = app.Simulation(prmtop.topology, system, integrator, platform,
                        platformProperties)
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)

    force_val_old = -spring_k*(keyInteraction_dist - dist_in)
    energy_val_old = u.Quantity(value=0, unit=u.kilocalorie/u.mole)

    f=open('duck_'+str(temperature).split()[0].replace('.0','K')+'_'+str(SMD_num)+'.dat','w')
    steps = int((dist_fin.value_in_unit(u.nanometer) / 0.000001 - dist_in.value_in_unit(u.nanometer) / 0.000001)) / move_force_step
    pull_distance = 0.000001 * move_force_step
    
    
    #write trajectory
    top = md.load_prmtop('system_solv.prmtop')
    atom_subset = top.select('not water')
    simulation.reporters.append(app.StateDataReporter("smd_"+str(temperature).split()[0].replace('.0','K')+"_"+str(SMD_num)+".csv", move_force_step, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True,
                                                      temperature=True, density=True, progress=True, totalSteps=move_force_step*steps, speed=True))
        
    simulation.reporters.append(HDF5Reporter("smd_"+str(temperature).split()[0].replace('.0','K')+"_"+str(SMD_num)+".h5", move_force_step*20, atomSubset=atom_subset))

    
    for i in range(steps):    
        state = simulation.context.getState(getPositions=True)
        pos_keyInt = state.getPositions()
        keyInteraction_pos = [pos_keyInt[keyInteraction[0]], pos_keyInt[keyInteraction[1]]]
    
        keyInteraction_dist = np.linalg.norm(keyInteraction_pos[0]-keyInteraction_pos[1])
        keyInteraction_vect = (keyInteraction_pos[1] - keyInteraction_pos[0]) / keyInteraction_dist
        keyInteraction_vect = u.Quantity(value=keyInteraction_vect, unit=u.nanometers)
        pullto = keyInteraction_pos[0] + (0.25 + float(i) * pull_distance) * keyInteraction_vect
    
        displace = pullto - pullto_old
    
        simulation.context.setParameter('displace_x', displace[0])
        simulation.context.setParameter('displace_y', displace[1])
        simulation.context.setParameter('displace_z', displace[2])
        if i == 0:
            distance = 0.0
        else:
            distance = pull_distance
        dist_spring =  (0.25 + float(i) * pull_distance) * u.nanometer
        force_val = -spring_k * (keyInteraction_dist - dist_spring)
        energy_val = energy_val_old + (distance * u.nanometer) * 0.5 * (force_val+force_val_old)
        force_val_old = force_val
        energy_val_old = energy_val
        if (i%int(save_step/move_force_step)) == 0:
            f.write(str(i)+' '+str(dist_spring)+' '+str(keyInteraction_dist)+' '+str(force_val)+' '+str(energy_val)+'\n')
        
        simulation.step(move_force_step)
        
    #f.write(str(i)+' '+str(dist_spring)+' '+str(keyInteraction_dist)+' '+str(force_val)+' '+str(energy_val)+'\n')
    f.close()
    

def get_Wqb_value(file_duck_dat, Plot):
    f = open(file_duck_dat,'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[1]), float(a[3]), float(a[5]), float(a[8])])
    f.close()
    data = np.array(data[1:])
    Work = data[:,3]
    #split it into segments of 200 points 
    num_segments = int(len(data)/200) 
    num_segments = int(len(data)/200) 
    #alayze each segment to see if minimum in the segment is the local minimum
    #local minimum is the point with the lowest value of 200 neighbouring points
    #first local minumum is miminum used later to duck analysis
    for segment in range(num_segments):
        #detecting minium inthe segment
        sub_data = data[segment * 200 : (segment + 1) * 200]
        sub_Work = sub_data[:,3] 
        index_local = np.argmin(sub_Work)
        #segment of 200 points arround detected minimum
        index_global = index_local + segment * 200
        if index_global > 100:
            sub2_data = data[index_global - 100 : index_global + 101]
        else:
            sub2_data = data[0 : index_global + 101]
        sub2_Work = sub2_data[:,3]
        index_local2 = np.argmin(sub2_Work)
        if index_global < 100:
            if index_local2 == index_global:
                
                Wqb_min_index = index_global
            break
        else:
            if index_local2 == 100:
                Wqb_min_index = index_global
                break
    
    Wqb_min = Work[Wqb_min_index]
    sub_max_data = data[Wqb_min_index:]
    sub_max_Work = sub_max_data[:,3]
    
    
    Wqb_max = max(sub_max_Work)
    
    Wqb_value = Wqb_max - Wqb_min
    if Plot == True:
        Dist_plot = data[:,0]*10
        Work_plot = Work - Wqb_min
    
        plt.plot(Dist_plot, Work_plot, 'b')
        plt.title(file_duck_dat)
        plt.xlabel('Dist [A]')
        plt.ylabel('Work [kcal/mol]')
        Wqb_max_index = np.argmax(sub_max_Work)
        Wqb_max_index_global = Wqb_max_index + Wqb_min_index
        position_max_x = Dist_plot[Wqb_max_index_global]
        position_max_y = Work_plot[Wqb_max_index_global]
        plt.annotate(str(Wqb_value), (position_max_x, position_max_y))
        plot_name = 'plot_'+file_duck_dat.replace('.dat','.png')
        plt.savefig(plot_name)
    return(Wqb_value)

###############################################################################    
    
#####################Sym parameters############################################

num_MD = 3 
dist_duck_rst_file = 'dist_duck.rst'   
spring_k = 100 * u.kilocalorie/(u.mole * u.angstrom * u.angstrom)
dist_in = 2.5 * u.angstrom
dist_fin = 5.0 * u.angstrom
save_step = 50
move_force_step = 50
Wqb_BLOCK_value = 6.0
Plot = False


###############################################################################    

######Conditional Import####
#the 

if Plot == True:
    import matplotlib.pyplot as plt 
    import seaborn as sns
    sns.set_style('ticks')

#Platform definition

platform = mm.Platform_getPlatformByName("CUDA")
platformProperties = {}
platformProperties['CudaPrecision'] = 'mixed'
platformProperties['CudaDeviceIndex'] = '0'

#Read files

prmtopName = 'system_solv.prmtop'
inpcrdName = 'system_solv.inpcrd'

prmtop = app.AmberPrmtopFile(prmtopName)
inpcrd = app.AmberInpcrdFile(inpcrdName)

#Minimization

print('Minimising...')

#Define system

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom)

Chunk_Heavy_Atoms = []

for atom_i in prmtop.topology.atoms():
    if atom_i.residue.name not in ('HOH', 'WAT', 'IP3', 'LIG', 'Cs+', 'K+', 'Rb+', 'Li+', 'Na+', 'IP', 'Cl-', 'IM', 'IB') and atom_i.name[0] != 'H':
        Chunk_Heavy_Atoms.append(atom_i.index)
        
keyInteraction = []

# See page 475 of http://ambermd.org/doc12/Amber17.pdf for format
f=open(dist_duck_rst_file, 'r')
for line in f:
    if line.strip().split()[0] == '&rst':
        key_atoms = line.strip().split()[1]
        key_atoms_1 = key_atoms.replace('iat=', '').replace(',', ' ').split()
        keyInteraction = [int(key_atoms_1[0])-1 , int(key_atoms_1[1])-1]
        
if keyInteraction == []:
    print('Did not select key interaction')
    
print('selected interaction: '+str(keyInteraction))

applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)

integrator = mm.VerletIntegrator(1*u.femtosecond)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(inpcrd.positions)

simulation.reporters.append(app.StateDataReporter("minimisation.csv", 1, step=True, potentialEnergy=True))

simulation.minimizeEnergy(maxIterations=1000)

#Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))

#Equlibration - heating
#new minimised positions, however using old restraints
print('Heating system under NVT...')

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds)
applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)
applyLigandChunkRestraint(system, 1.0, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)
velocities = ()

for temperature in (150.0, 200.0, 250.0, 300.0):
    
    print('...to '+str(temperature).replace('.0', '')+' K...')

    integrator = mm.LangevinIntegrator(temperature*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)
    simulation.context.setPositions(positions)   
    if temperature != 150.0:
        simulation.context.setVelocities(velocities)

    simulation.reporters.append(app.StateDataReporter('heating'+str(temperature).replace('.0', '')+'.csv', 2000, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, density=True, progress=True, totalSteps=50000, speed=True))
    simulation.step(50000)
    #get positions for the next step of symulation
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions()
    velocities = state.getVelocities()
    
app.PDBFile.writeFile(simulation.topology, positions, open('heating_final.pdb', 'w'))

#Equlibration - density
'''
In AMBER version were using Berendsen barostat
The Berendsen barostat is not impemented in OpenMM thus not applied here.
Mone Carlo is the only method here.

'''
    
print('Density correction under NPT...')  # Now that we are warm, let's allow the density to correct for 70 ps

integrator = mm.LangevinIntegrator(temperature*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)

barostat = mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin, 25)
system.addForce(barostat)

simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)
simulation.reporters.append(app.StateDataReporter("density.csv", 2000, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True,
                                                  temperature=True, density=True, progress=True, totalSteps=500000, speed=True))

simulation.step(500000)

# Save the positions and velocities
state = simulation.context.getState(getPositions=True, getVelocities=True)
positions = state.getPositions()
velocities = state.getVelocities()
app.PDBFile.writeFile(simulation.topology, positions, open('density_final.pdb', 'w'))


#New system without ligand chunk restraints 

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds)
applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)

#SMD_0
print('SMD_0')
SMD_num = 0
temperature = 300.0 * u.kelvin
SMD(system, prmtop, platform, platformProperties, temperature, positions, velocities, keyInteraction, spring_k, dist_in, dist_fin, SMD_num, save_step, move_force_step)

#BLOCK
Wqb_value = get_Wqb_value('duck_'+str(temperature).split()[0].replace('.0','K')+'_'+str(SMD_num)+'.dat', Plot)
print('Wqb value '+str(Wqb_value))
if Wqb_value > Wqb_BLOCK_value:

    #SMD_325K_0
    print('SMD_325K_0')
    temperature = 325.0 * u.kelvin
    SMD(system, prmtop, platform, platformProperties, temperature, positions, velocities, keyInteraction, spring_k, dist_in, dist_fin, SMD_num, save_step, move_force_step)

    #BLOCK
    Wqb_value = get_Wqb_value('duck_'+str(temperature).split()[0].replace('.0','K')+'_'+str(SMD_num)+'.dat', Plot)
    print('Wqb value '+str(Wqb_value))
    if Wqb_value > Wqb_BLOCK_value:

        #Loop for N
        for MD in range(num_MD):

            #MD_N
            MD_num = MD+1
            print('MD_'+str(MD_num))
    
            system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds)
            applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)
            applyLigandChunkRestraint(system, 1.0, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)

            integrator = mm.LangevinIntegrator(300.0*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)
            simulation = app.Simulation(prmtop.topology, system, integrator, platform,
                                        platformProperties)
            simulation.context.setPositions(positions)
            simulation.context.setVelocities(velocities)
    
            simulation.reporters.append(app.StateDataReporter("md_"+str(MD_num)+".csv", 2000, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True,
                                                              temperature=True, density=True, progress=True, totalSteps=500000, speed=True))


            simulation.step(250000)
            
            state = simulation.context.getState(getPositions=True, getVelocities=True)
            positions = state.getPositions()
            velocities = state.getVelocities()
            app.PDBFile.writeFile(simulation.topology, positions, open('md_'+str(MD_num)+'.pdb', 'w'))
    
            #SMD
    
            system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds)
            applyHarmonicPositionalRestraints(system, 1.0, inpcrd, Chunk_Heavy_Atoms)

    
            #SMD_N
            SMD_num = MD+1
            print('SMD_300K_'+str(SMD_num))
            temperature = 300.0 * u.kelvin
            SMD(system, prmtop, platform, platformProperties, temperature, positions, velocities, keyInteraction, spring_k, dist_in, dist_fin, SMD_num, save_step, move_force_step)
            Wqb_value1 = get_Wqb_value('duck_'+str(temperature).split()[0].replace('.0','K')+'_'+str(SMD_num)+'.dat', Plot)
            print('Wqb value '+str(Wqb_value1))
            
            #SMD_325K_N
            print('SMD_325K_'+str(SMD_num))
            temperature = 325.0 * u.kelvin
            SMD(system, prmtop, platform, platformProperties, temperature, positions, velocities, keyInteraction, spring_k, dist_in, dist_fin, SMD_num, save_step, move_force_step)
            Wqb_value2 = get_Wqb_value('duck_'+str(temperature).split()[0].replace('.0','K')+'_'+str(SMD_num)+'.dat', Plot)
            print('Wqb value '+str(Wqb_value2))
            
            #BLOCK
            if Wqb_value1 > Wqb_BLOCK_value and Wqb_value2 > Wqb_BLOCK_value:
                continue
            else:
                print('Wqb value lower than '+str(Wqb_BLOCK_value))
                break
            
    else:
        print('Wqb value lower than '+str(Wqb_BLOCK_value))
else:
    print('Wqb value lower than '+str(Wqb_BLOCK_value))

print('DONE')
