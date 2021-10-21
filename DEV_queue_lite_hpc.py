

import os 
import time
import shutil
import random
from progress.bar import Bar

import DEV_graphene_maker_hpc as gmak


###### USER INPUTS, see README file. #######
dir_path = os.path.dirname(os.path.realpath(__file__))
lammps_files_path = "%s/LAMMPS_files"%dir_path
input_file_name = 'in.multi_bombard'
data_file_name = 'data.graphite_sheet'
results_dir_name = 'results'
#energies = [7, 10, 13, 16, 20, 25, 30, 40, 50, 65, 80, 100] #7ev minimum
energies = [100]
#repeats = 50
pre_bomb_run_val = '3000' #this is changed if test == True, so must be changed here
                            #rather than the input file.
bomb_run_val = "500" #this is automatically increase for slow particles and must be
                    #changed here rather than the input file.
number_of_particles = '375' #bombarding
test = False
hpc = True    
if test == True:
    number_of_particles = '5'
    #repeats = 3
############################################

full_start = time.time()

in_file = gmak.file_proc("%s/%s"%(lammps_files_path, input_file_name))
#velocities = ['200'] #Initial Low energy particle used to find surface
surface_finder = True                

all_files = []
results_list = []
paths = []
progress_markers_paths = []
replicate = None
#for velocity in velocities: #loops through velocities correspoinding to inputted energies.
                            #These are calculated and appended later.
repeat_count = 0
directory_names = []

#while repeat_count < repeats: #loops through repeat simulations for the same energy
    
new = ''
for i in in_file: #goes through the input file line by line both reading and editing
    index = in_file.index(i)
    i = i.split()
    seperator = ' '

    try:
        if i[0] == '#number_of_sheets':
            graphite_sheets = float(i[1])

        if i[0] == '#atom_type': #either 'h', 'd' or 't'
            atom_type = i[1]

        if i[0] == "variable":
            i[-1] = str(number_of_particles)
            in_file[index] = seperator.join(i)     

        if i[0] == 'read_data':

            if i[1] == 'data.diamond':
                i[1] = "%s/data.diamond"%lammps_files_path
            if i[1] == 'data.graphite_sheet':
                i[1] = "%s/data.graphite_sheet"%lammps_files_path
        
            in_file[index] = seperator.join(i)         

        if i[0] == 'pair_coeff':

            i[3] = "%s/CH.rebo"%lammps_files_path
            in_file[index] = seperator.join(i)

        if i[0] == "set": #setting velocity line
            #if len(velocities) == 1: #make
            if atom_type == 'h':
                atom_mass = 1.0079
            if atom_type == 'd':
                atom_mass = 2.0014
            if atom_type == 't':
                atom_mass = 3.0160
            if atom_type == 'ten':
                atom_mass = 10
            if atom_type != 'h' and atom_type != 't' and atom_type !='d':
                atom_mass = float(atom_type)

            energy = energies[0]
            unit_conv = 1e-2*(1.602e-19/1.661e-27)**0.5
            velocity = ((2*energy/atom_mass)**0.5)*unit_conv 

            i[-1] = str(velocity)
            in_file[index] = seperator.join(i)
        
        
        if i[0] == "replicate":
            replicate = i[1:]
            """
            if float(velocity) > 1100:
                i[-1] = '10'
            else:
                i[-1] = '8'
            """
            in_file[index] = seperator.join(i)
        

        if i[0] == "fix" and len(i[-1]) == 5: #changing random seed for each repeat
            rand = random.randint(10000,99999)
            i[-1] = str(rand)
            in_file[index] = seperator.join(i)

        if i[0] == "velocity" and i[1] == 'all' and i[2] == 'create': #changing random seed for each repeat
            rand = random.randint(10000,99999)
            i[4] = str(rand)
            in_file[index] = seperator.join(i)

        if i[0] == 'region' and i[1] == 'box': #creating region for bombardment atom to be created
            central = True                      #this is adjusted for the particular simulation
            diamond_size = 3.567
            graphene_thickness = 3.35

            #if replicate == None: required if using one data file for diamond and no replicate
             #   replicate = [8,8,10]

            if central == True:
                xlo, ylo = [float(replicate[i])*diamond_size*(-1) for i in range(0,2)]
                xhi, yhi = [(float(replicate[i]) + 1)*diamond_size for i in range(0,2)]

            else:
                xlo = 0
                xhi = float(replicate[0])*diamond_size
                ylo = 0
                yhi = float(replicate[1])*diamond_size

            zhi = -graphene_thickness*graphite_sheets - 25
            zlo = zhi - 30
            
            i = seperator.join(i[:3]) + " %s %s %s %s %s %s"%(xlo, xhi, ylo ,yhi,zlo,zhi)

            in_file[index] = i

        if i[0] == 'run' and i[-1] == '#':
            if float(velocity) < 300:
                i[1] = '1000'
            else:
                i[1] = bomb_run_val

            in_file[index] = seperator.join(i)

        if i[0] == 'run' and i[-1] != '#':
            if test == True:
                i[1] = '3000'
            else:
                i[1] = pre_bomb_run_val

            in_file[index] = seperator.join(i)

    except IndexError:                                    
        pass

    new += "%s\n"%in_file[index]

with open("%s/%s"%(lammps_files_path,input_file_name), 'w') as fp: #rewriting edited input file
    fp.write(str(new)) 

if atom_type == 'h':
    atom_mass = 1.0079
if atom_type == 'd':
    atom_mass = 2.0014
if atom_type == 't':
    atom_mass = 3.0160

#if len(velocities) == 1: #adds corresponding velocities from energy inputs
#    unit_conv = 1e-2*(1.602e-19/1.661e-27)**0.5
#    new_velocities = [((2*energy/atom_mass)**0.5)*unit_conv for energy in energies] 
#    velocities += [str(round(velocity,3)) for velocity in new_velocities]

#writes graphene data sheet for the specific simulation
#returns dictionary containing the number of bulk atoms in each region
#used for lattice displacement calculations

if surface_finder == True: #only run onces
    bulk_atoms_dict = gmak.main(path = lammps_files_path, data_file_name = 'data.graphene', input_file_name = input_file_name, 
                                new_data_file_name = data_file_name, atom_mass = atom_mass)

count = 1
while True: #creating new directory 
    if surface_finder == True:
        velocity = f"{energy}eV"
    try:
        new_path = "%s/%s/%s_%sg_%s_%s_%s"%(dir_path, results_dir_name, atom_type, int(graphite_sheets), 
                                            velocity, number_of_particles, count)
        os.mkdir(new_path)
        paths.append(new_path)
        if count == 1:
            progress_markers_paths.append(new_path)
        break
    except FileExistsError:
        count += 1  

shutil.copyfile("%s/data.graphite_sheet"%lammps_files_path, "%s/data.graphite_sheet"%new_path)
shutil.copyfile("%s/%s"%(lammps_files_path, input_file_name), "%s/%s"%(new_path, input_file_name))

if hpc == True:         
    shutil.copyfile("srun","%s/srun"%new_path) #copying srun avoids rewritting
    os.chdir(new_path)
    os.system("sbatch %s/srun"%new_path)  #submitted job to DIRAC

lammps_files_path = "%s/LAMMPS_files"%dir_path

if hpc == False: #if run non-parallel, although standard package is suggested
    os.chdir(new_path)    
    os.system("mpiexec -n 2 lmp_serial -in %s/%s"%(lammps_files_path, input_file_name))

#repeat_count += 1

#if surface_finder == True: #No repeats for initial surface probe run
 #   surface_finder = False
  #  break


with open("%s/%s/file_paths.txt"%(dir_path, results_dir_name), 'w') as fp: #rewriting edited input file
    fp.write(str(paths)) 

progress_markers_paths.append(new_path)
start = time.time()
with Bar('Processing', max=len(progress_markers_paths)) as bar: #progress bar

    for path in progress_markers_paths:
        counter = 0
        while True: #loop searches for log.lammps file to update progress
            time.sleep(5)
            counter += 5
            try:
                open("%s/log.lammps"%path)
                break
            except FileNotFoundError:
                pass

            if counter > 600:
                print("Waited for 10minutes. No log file yet. \nNot looking in %s anymore."%path)
                break
        bar.next()

timer = time.time() - start
time.strftime('%H:%M:%S', time.gmtime(timer))
print("Final log file recieved.")
print("This took %s"%time.strftime('%H:%M:%S', time.gmtime(timer)))


