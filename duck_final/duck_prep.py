# function for getting the data
def get_value(data, value):
    for line in data:
        if line[:len(value)] == value:
            val_ret = line.replace(' ','').replace('=',' ').split()[1]
            break
    return(val_ret)

def mod_list(open_file, line_temp, line_new):
    modif_list = []
    for line in open_file:
        if line == line_temp:
            modif_list.append(line_new)
        else:
            modif_list.append(line)
    return(modif_list)

# read setup.duck file to get parameters
f = open('setup.duck','r')
setup_duck = []
for line in f :
    setup_duck.append(line.strip())
f.close()

# Get variables

prmtop = get_value(setup_duck,'prmtop')
inpcrd = get_value(setup_duck,'inpcrd')
key_int = get_value(setup_duck,'key_int').replace(',',' ').split()
MD_traj_frames = get_value(setup_duck,'MD_traj_frames')
SMD_traj_frames = get_value(setup_duck,'SMD_traj_frames') 
max_SMD_runs = get_value(setup_duck,'max_SMD_runs')
Wqb_thresh = get_value(setup_duck,'Wqb_thresh')
MD_len = get_value(setup_duck,'MD_len')
SMD_dist_start = get_value(setup_duck,'SMD_dist_start')
SMD_dist_end = get_value(setup_duck,'SMD_dist_end')
k_str = get_value(setup_duck,'k_str')

# Set simulation files
# equil.py

equil_file = []
f = open('temp_equil.py','r')
for line in f:
    equil_file.append(line)
f.close()

equil_file = mod_list(equil_file, "#!#keyInteraction_ind_mol2 = [XXX, YYY]\n", "keyInteraction_ind_mol2 = [" + key_int[0] + ", " + key_int[1] + "]\n")
equil_file = mod_list(equil_file, "#!#prmtopName = 'PRMTOP'\n", "prmtopName = '"+prmtop+"'\n")
equil_file = mod_list(equil_file, "#!#inpcrdName = 'INPCRD'\n", "inpcrdName = '"+inpcrd+"'\n")

f = open('equil.py', 'w')
for line in equil_file:
    f.write(line)
f.close()

# md.py

md_file = []
f = open('temp_md.py','r')
for line in f:
    md_file.append(line)
f.close()

md_file = mod_list(md_file, "#!#prmtopName = 'PRMTOP'\n", "prmtopName = '"+prmtop+"'\n")
md_file = mod_list(md_file, "#!#inpcrdName = 'INPCRD'\n", "inpcrdName = '"+inpcrd+"'\n")
md_file = mod_list(md_file, "#!#keyInteraction_ind_mol2 = [XXX, YYY]\n", "keyInteraction_ind_mol2 = [" + key_int[0] + ", " + key_int[1] + "]\n")
md_file = mod_list(md_file,'#!#MD_len = X * u.nanosecond\n','MD_len = '+MD_len+' * u.nanosecond\n')
if MD_traj_frames == '0':
    md_file = mod_list(md_file, '#!#from mdtraj.reporters import HDF5Reporter\n','\n')
    md_file = mod_list(md_file, '#!#traj_out_file = sys.argv[5]\n','\n')
    md_file = mod_list(md_file, '#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))\n','\n')
else:
    md_file = mod_list(md_file, '#!#from mdtraj.reporters import HDF5Reporter\n','from mdtraj.reporters import HDF5Reporter\n')
    md_file = mod_list(md_file, '#!#traj_out_file = sys.argv[5]\n','traj_out_file = sys.argv[5]\n')
    md_file = mod_list(md_file, '#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))\n','simulation.reporters.append(HDF5Reporter(traj_out_file, '+MD_traj_frames+'))\n')
   
f = open('md.py', 'w')
for line in md_file:
    f.write(line)
f.close()


# smd.py

smd_file = []
f = open('temp_smd.py','r')
for line in f:
    smd_file.append(line)
f.close()

smd_file = mod_list(smd_file, "#!#prmtopName = 'PRMTOP'\n", "prmtopName = '"+prmtop+"'\n")
smd_file = mod_list(smd_file, "#!#inpcrdName = 'INPCRD'\n", "inpcrdName = '"+inpcrd+"'\n")
smd_file = mod_list(smd_file, "#!#keyInteraction_ind_mol2 = [XXX, YYY]\n", "keyInteraction_ind_mol2 = [" + key_int[0] + ", " + key_int[1] + "]\n")
smd_file = mod_list(smd_file, '#!#MD_len = X * u.nanosecond\n','MD_len = '+MD_len+' * u.nanosecond\n')
smd_file = mod_list(smd_file, "#!#spring_k = XX * u.kilocalorie/(u.mole * u.angstrom * u.angstrom)\n", "spring_k = "+k_str+" * u.kilocalorie/(u.mole * u.angstrom * u.angstrom)\n")
smd_file = mod_list(smd_file, "#!#dist_in = XXX * u.angstrom # in angstrom\n", "dist_in = "+SMD_dist_start+" * u.angstrom # in angstrom\n")
smd_file = mod_list(smd_file, '#!#dist_fin = XXX * u.angstrom # in angstrom\n','dist_fin = '+SMD_dist_end+' * u.angstrom # in angstrom\n')


if SMD_traj_frames == '0':
    smd_file = mod_list(smd_file, '#!#from mdtraj.reporters import HDF5Reporter\n','\n')
    smd_file = mod_list(smd_file, '#!#traj_out_file = sys.argv[6]\n','\n')
    smd_file = mod_list(smd_file, '#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))\n','\n')
else:
    smd_file = mod_list(smd_file, '#!#from mdtraj.reporters import HDF5Reporter\n','from mdtraj.reporters import HDF5Reporter\n')
    smd_file = mod_list(smd_file, '#!#traj_out_file = sys.argv[6]\n','traj_out_file = sys.argv[6]\n')
    smd_file = mod_list(smd_file, '#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))\n','simulation.reporters.append(HDF5Reporter(traj_out_file, '+SMD_traj_frames+'))\n')
   
f = open('smd.py', 'w')
for line in smd_file:
    f.write(line)
f.close()

# job.q

f = open('job.q','w')
f.write('#!/bin/bash\n#SBATCH -J OMM-7\n#SBATCH --gres=gpu:1\n#SBATCH --constraint=m2090\n\n')
f.write('module purge\nmodule add lapack/gcc/64/3.5.0\nmodule add openblas/dynamic/0.2.14\nmodule add numpy/1.8.1\nmodule add matplotlib/1.3.1\nmodule add scipy/0.16.1\nmodule add hdf5_18/1.8.14\nmodule add pytables/3.3.0\nmodule add pandas/0.19.2\nmodule add mdtraj/1.8.0\nmodule add cuda75/toolkit\nmodule add openmm/7.0.1\n\n')
f.write('#equil\n')
f.write('python equil.py > equil.log\n\n')
f.write('#SMD 300K 0\n')
if SMD_traj_frames == '0':
    f.write('python smd.py 300 equil.chk smd_300_0.csv smd_300_0.dat smd_300_0.pdb > smd_300_0.log\n')
else:
    f.write('python smd.py 300 equil.chk smd_300_0.csv smd_300_0.dat smd_300_0.pdb smd_300_0.nc > smd_300_0.log\n')
f.write('\n#Wqb BLOCK\n')
f.write('Wqb=$(python getWqbValue.py)\n')
f.write("flag=$(echo $Wqb | awk '{if($1>"+Wqb_thresh+") print 1; else print 0}')\n")
f.write('if [ $flag -eq 0 ]\n')
f.write('then\n')
f.write('exit\n')
f.write('fi\n\n')
f.write('#SMD 325K 0\n')
if SMD_traj_frames == '0':
    f.write('python smd.py 325 equil.chk smd_325_0.csv smd_325_0.dat smd_325_0.pdb > smd_325_0.log\n')
else:
    f.write('python smd.py 325 equil.chk smd_325_0.csv smd_325_0.dat smd_325_0.pdb smd_325_0.nc > smd_325_0.log\n')

f.write('\n#Wqb BLOCK\n')
f.write('Wqb=$(python getWqbValue.py)\n')
f.write("flag=$(echo $Wqb | awk '{if($1>"+Wqb_thresh+") print 1; else print 0}')\n")
f.write('if [ $flag -eq 0 ]\n')
f.write('then\n')
f.write('exit\n')
f.write('fi\n\n')

for num in range(max_SMD_runs):
    NUM=num+1
    f.write('#md'+str(NUM)+'\n')
    if num == 0:
        if MD_traj_frames == '0':
            f.write('python md.py equil.chk md'+str(NUM)+'.chk md'+str(NUM)+'.csv md'+str(NUM)+'.pdb > md'+str(NUM)+'.log\n\n')
        else:
            f.write('python md.py equil.chk md'+str(NUM)+'.chk md'+str(NUM)+'.csv md'+str(NUM)+'.pdb md'+str(NUM)+'.nc > md'+str(NUM)+'.log\n\n')
    else:
        if MD_traj_frames == '0':
            f.write('python md.py md'+str(num)+'.chk md'+str(NUM)+'.chk md'+str(NUM)+'.csv md'+str(NUM)+'.pdb > md'+str(NUM)+'.log\n\n')
        else:
            f.write('python md.py md'+str(num)+'.chk md'+str(NUM)+'.chk md'+str(NUM)+'.csv md'+str(NUM)+'.pdb md'+str(NUM)+'.nc > md'+str(NUM)+'.log\n\n')
    
    f.write('#SMD 300K '+str(NUM)+'\n')
    if SMD_traj_frames == '0':
        f.write('python smd.py 300 md'+str(NUM)+'.chk smd_300_'+str(NUM)+'.csv smd_300_'+str(NUM)+'.dat smd_300_'+str(NUM)+'.pdb > smd_300_'+str(NUM)+'.log\n\n')
    else:
        f.write('python smd.py 300 md'+str(NUM)+'.chk smd_300_'+str(NUM)+'.csv smd_300_'+str(NUM)+'.dat smd_300_'+str(NUM)+'.pdb smd_300_'+str(NUM)+'.nc > smd_300_'+str(NUM)+'.log\n\n')
    f.write('#Wqb BLOCK\n')
    f.write('Wqb=$(python getWqbValue.py)\n')
    f.write("flag=$(echo $Wqb | awk '{if($1>"+Wqb_thresh+") print 1; else print 0}')\n")
    f.write('if [ $flag -eq 0 ]\n')
    f.write('then\n')
    f.write('exit\n')
    f.write('fi\n\n')
    f.write('#SMD 325K '+str(NUM)+'\n')
    if SMD_traj_frames == '0':
        f.write('python smd.py 325 md'+str(NUM)+'.chk smd_325_'+str(NUM)+'.csv smd_325_'+str(NUM)+'.dat smd_325_'+str(NUM)+'.pdb > smd_325_'+str(NUM)+'.log\n\n')
    else:
        f.write('python smd.py 325 md'+str(NUM)+'.chk smd_325_'+str(NUM)+'.csv smd_325_'+str(NUM)+'.dat smd_325_'+str(NUM)+'.pdb smd_325_'+str(NUM)+'.nc > smd_325_'+str(NUM)+'.log\n\n')
    
    f.write('\n#Wqb BLOCK\n')
    f.write('Wqb=$(python getWqbValue.py)\n')
    f.write("flag=$(echo $Wqb | awk '{if($1>"+Wqb_thresh+") print 1; else print 0}')\n")
    f.write('if [ $flag -eq 0 ]\n')
    f.write('then\n')
    f.write('exit\n')
    f.write('fi\n\n')

f.write('python getWqbValuePLOT.py > wqb_value.txt\n')

f.close()

