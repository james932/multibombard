

import os
import shutil
import sys
import matplotlib.pyplot as plt

import DEV_log_interpretter_lite_hpc as lint
import DEV_atom_track_multi as tracker
import DEV_graphene_maker_hpc as gmak
import DEV_result_plotter_hpc as plot



###TODO#####

    # Check box sizes for geting further than 215th atom
    # Bulk atom checks (displacement) - put limit on so 40angstrom movements are counted
    # Get .xyz file for last prebombardment step 
    # Make get_frame.py that extracts frame form all.xyz

############
'''
This script is an alteration on the original analyse.py. The main difference is in the
structure. The original scrip was designed to be run directly after a LAMMPS sim, whereas
this script is designed to be run on larger collection of simualtions. This seperates
queue_hpc.oy and this script. So all jobs can be sent off to DIRAC and when completed,
this script can be run to analyse. Removing the issue of the varying amount of time
taken for a job to return from DIRAC.
'''

class Vel_Dir:
    """
    A directory object. This refers to the directory containing all of the 
    """
    def __init__(self, velocity, sim_dirs):
        self.velocity = velocity
        self.sim_dirs = sim_dirs
        self.surface_finder = False
    
    def update(self):
        if len(self.sim_dirs) == 1:
            self.surface_finder = True

        self.energy = self.sim_dirs[0].energy
        self.path = self.sim_dirs[0].vel_dir_path

class Sim_Directory:

    def __init__(self, path, folder_name, repeat_folder_name):

        self.folder_name = folder_name
        self.repeat_folder_name = repeat_folder_name
        self.path = path
        self.full_path = "%s/%s/%s"%(path,folder_name,repeat_folder_name)
        self.vel_dir_path = "%s/%s"%(path,folder_name)
        self.surface_finder = False

        split_name = folder_name.split("_") 

        self.velocity = split_name[-1]
        self.atom_type = split_name[0]

        if self.velocity == 'probe':
            self.surface_finder = True
            self.energy = 5

        else:
            self.velocity = float(self.velocity)
            self.energy = self.energy_calc()

    def energy_calc(self):

        if self.atom_type == 'h':
            atom_mass = 1.0079
        if self.atom_type == 'd':
            atom_mass = 2.0014
        if self.atom_type == 't':
            atom_mass = 3.0160
        if self.atom_type == 'ten':
            atom_mass = 10
        if self.atom_type != 'h' and self.atom_type != 't' and self.atom_type !='d':
            atom_mass = float(self.atom_type)

        unit_conv = 1.661e-27*1e4/1.602e-19
        energy = (0.5*atom_mass*self.velocity**2)*unit_conv

        return energy

def find_paths():
    
    paths_txt = open("file_paths.txt", 'r')
    paths_txt = paths_txt.read()
    paths_txt = paths_txt[2:-2]
    paths_list = paths_txt.split("', '")

    return paths_list

def fail_checker(sim_dirs):

    failed_sims = []
    for sim_dir in sim_dirs:
        try:
            paths_txt = open("%s/FAIL.txt"%sim_dir.full_path, 'r')
            setattr(sim_dir, 'success', False)
            failed_sims.append(sim_dir.repeat_folder_name)

        except FileNotFoundError:
            success = True
            setattr(sim_dir, 'success', True)

    return failed_sims

    

def sort(paths_list):

    
    paths_to_sim = []

    for result_dir in paths_list:

        lint.main(path = result_dir, in_file= "%s/in.multi_bombard"%result_dir, 
                data_file = "%s/data.graphite_sheet"%result_dir, reduce_size= True)

        path = result_dir.split('/')
        dash  = '_'
        folder_name = path[-1]
        repeat_folder_name = folder_name
        folder_name = folder_name.split('_')
        velocity_folder_name = folder_name[:-1]
        velocity_folder_name = dash.join(velocity_folder_name)

        slash = '/'
        path = slash.join(path[:-1])

        while True:
            try:
                os.mkdir("%s/%s"%(path,velocity_folder_name))

            except FileExistsError:
                shutil.move(result_dir, "%s/%s"%(path,velocity_folder_name))
                paths_to_sim.append(Sim_Directory(path,velocity_folder_name, repeat_folder_name))
                break
    
    return paths_to_sim


    
def main():

    current_dir = os.path.dirname(os.path.realpath(__file__))
    path_end = sys.argv[1]

    path = current_dir + '/hpc_results/' + path_end

    path_end_split = path_end.split('_')

    for i in path_end_split:
        try:
            if i[-2:] == "eV":
                energy = float(i[:-2])
        except IndexError:
            pass


    lint.main(path = path, in_file= '%s/in.multi_bombard'%path, data_file = '%s/data.graphite_sheet'%path, reduce_size = True)

    bulk_atoms_dict, graphite_sheets, atom_type, diamond_type, replicate, no_bombarding_atoms = gmak.main(path,
                                                                                            data_file_name=None, 
                                                                                            input_file_name="in.multi_bombard", 
                                                                                            new_data_file_name= None)
    print(f"anayls atom_type: {atom_type} ")

    results = tracker.main(all_xyz_file_path = "%s/all.xyz"%path, no_of_sheets = graphite_sheets, atom_type = atom_type, 
                            diamond_type = diamond_type, energy = energy, replicate = replicate, no_bombarding_atoms = no_bombarding_atoms)


    count = 1
    while True:
        try:
            new_path = "%s/%s_%s_layers_%s"%(path, atom_type, int(graphite_sheets), count)
            os.mkdir(new_path)
            break
        except FileExistsError:
            count += 1  



    with open("%s/results.txt"%new_path, 'w') as fp:
        fp.write(results.present())

    plt.hist(results.surface_depth_profile, density = True, stacked = True)
    plt.xlabel("Penetration / Å ")
    plt.ylabel("Fraction of Retained Particles")
    #plt.vlines(x = 0,ymin = 0, ymax= 0.5, linestyles='dashed')
    plt.savefig("%s/depth_profile.png"%new_path)
    plt.show()

main()


