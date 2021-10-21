

import numpy as np
from collections import namedtuple

'''
- Need to edit quite a lot for multi bomdard
- Want number of bobmarded atoms vs time
- Added some time after the last bombardment, then can just count the number of
  2 atoms at the end.
- Average depth
- Max depth
- Some comparison between before and after carbon atoms - will need some work

'''


class Particle:
    """Object that calculates and stores path data and results."""

    def __init__(self, repeat_no, path_str, no_of_sheets, surface = 0, bulk = False):
        """path_str is a list of positions in string format pulled directly
        from all.xyz.
        bulk represents wether the atom is part of the bulk material
        (limiting calculations in this case) or if it is a bombarding particle. 
        surface defaults to 0 but is taken as an argument once the surface height
        has been determined by the initial surface probe run."""

        self.repeat_no = repeat_no
        self.bulk = bulk

        self.xs = [float(pos[0]) for pos in path_str]
        self.ys = [float(pos[1]) for pos in path_str]
        self.zs = [float(pos[2]) for pos in path_str]

    
        if no_of_sheets == 0:
            surface = 0

        self.max_penetration = max(self.zs) - surface

        if max(self.zs) > 0:
            self.diamond_max_penetration = max(self.zs)
        else:
            self.diamond_max_penetration = 0

        self.path_float = [[float(pos[0]),float(pos[1]),float(pos[2])] for pos in path_str]

        self.path_string = ''
        for i in path_str:
            self.path_string += str(i) + '\n'

        self.bounce = False
        self.diamond_pen = False
        self.rebound(surface)
       
        self.string = self.present(path_str)


    def present(self, path_str):
        """Produces string of information for each bombardment particle."""

        repeat = '\nRepeat Number: ' + str(self.repeat_no)
        penetration = "\nMax Penetration: " + str(self.max_penetration)
        path = "\nParticle Path:\n" + str(self.path_string) + '\n'
        
        if self.bounce != None:
            bounce_str = "\nThe particle was reflected from the %s."%self.bounce
        else:
            bounce_str = "\nThe particle was not reflected"

        if self.pen != None:
            pen_str = "\nThe particle penetrated the %s."%self.pen
        else:
            pen_str = "\nThe particle did not penetrate the material."

        if self.finish != None:
            finish_str = "\nThe particles final location: %s."%self.finish
        else:
            finish_str = "\nThe particles final location could not be determined.."


        return repeat + penetration + bounce_str + pen_str + finish_str + path

    def rebound(self, surface):
        """Calculates and attributes penetration regions, reflections (bounce) 
        and final locations. """

        self.pen = None
        self.bounce = None
        self.finish = None

        if self.bulk == False:

            grad = np.diff(self.zs)
            for i in range(0,len(grad)):

                if self.zs[i] > (surface + 0.2):
                    self.pen = "graphene"

                if self.zs[i] > 0.2:
                    self.pen = "graphene and diamond"

                if grad[i] < 0 and self.zs[i] < surface:
                    self.bounce = "graphene"

                if grad[i] < 0 and self.zs[i] > surface and self.zs[i] < 0:
                    self.bounce = "diamond"

        if self.zs[-1] < 0 and self.zs[-1] >= surface:
            self.finish = "graphene"

        if self.zs[-1] < surface:
            self.finish = "outside"

        if self.zs[-1] >= 0:
            self.finish = "diamond"

            


def file_proc(file):
    """Short function to split text files into a list of lines"""

    opened_file = open(file, 'r')
    opened_file = opened_file.read()
    opened_file = opened_file.split("\n")
    return opened_file


def main(path, directory_names, no_of_sheets, atom_type, energy, diamond_type, surface = 0, bulk_atoms_dict = None):

    item = namedtuple('item', ['name', 'val'])
    atom = namedtuple('atom', ['id', 'type', 'pos'])
    
    if type(directory_names[0]) != str:
        directory_names = [sim_file.full_path for sim_file in directory_names if sim_file.success == True]

    particles  = [] #bombarding
    bulk_particles = []
    bombardment_particles = []
    items = []

    #used for graph titles.
    if atom_type == 'h': 
        atom_type = "Hydrogen"
    if atom_type == 'd':
        atom_type = "Deuterium"
    if atom_type == 't':
        atom_type = "Tritium"
        
    for directory in directory_names: #loops through all repeats of an energy

        all_xyz = '%s/all.xyz'%(directory)

        try:
            all_xyz = file_proc(all_xyz)
        except FileNotFoundError:
            print("\nERROR: %s, File Not found."%all_xyz)
            break

        bombardment_path = []
        bulk_particle = []

        time_split = all_xyz.split('\n\n')

        for time_step in time_split:
            time_step = time_step.split('\n')
            items = []

            for i in time_step:
                
                
        #looping through lines in xyz file
                i = i.split()
                seperator = ' '

                try:
                    if i[0] == 'ITEM:':
                        item_bool = True
                        item_title = i[1:]
                    
                    if item_bool == True:
                        item = item(name = item_title, val = i)
                        items.append(item)
                        item_bool = False

                    #if i[0] == 'Atoms.':
                     #   final_file_line_index = all_xyz.index(seperator.join(i))
                        
                    if i[1] == '1':
                        bulk_particles.append(atom(id = i[0], type = i[1], pos = i[2:]))
                    
                    if i[1] == '2': #tracking bombardment atom
                        bombardment_particles.append(atom(id = i[0], type = i[1], pos = i[2:], items = items))
                except IndexError:
                    pass

            for i in range(final_file_line_index+1, len(all_xyz)): #only for last time step

                line = all_xyz[i]
                line = line.split()

                try:
                    if line[0] == "1": #bulk atoms
                        bulk_particle.append(Particle(repeat_no = int(directory_names.index(directory)+1), path_str= [line[1:]],
                                                        no_of_sheets = no_of_sheets, surface = surface, bulk = True))
                except IndexError:
                    pass

        #bombarding particles, each particle is a repeat.
        particles.append(Particle(repeat_no = int(directory_names.index(directory)+1), path_str = bombardment_path, 
                                    no_of_sheets = no_of_sheets, surface = surface))
        
        #number of bulk particles in each region
        bulk_atoms_diamond = len([particle.finish for particle in bulk_particle if particle.finish == 'diamond'])
        bulk_atoms_graphene = len([particle.finish for particle in bulk_particle if particle.finish == 'graphene'])
        bulk_atoms_outside = len([particle.finish for particle in bulk_particle if particle.finish == 'outside'])

        bulk_particles.append(dict(bulk_atoms_diamond = bulk_atoms_diamond, bulk_atoms_graphene = bulk_atoms_graphene,
                                        bulk_atoms_outside = bulk_atoms_outside))
        

    #averaging over repeats for one energy
    #bulk_atoms_dict is produced by graphene_maker.py and gives number of bulk atoms 
    #in each region initially.
        
    bulk_diamond_vals = [repeat['bulk_atoms_diamond'] for repeat in bulk_particles]
    bulk_diamond_avg_err = [np.average(bulk_diamond_vals), np.std(bulk_diamond_vals)]
    bulk_diamond_lost = bulk_atoms_dict["diamond"] - bulk_diamond_avg_err[0]

    bulk_graphene_vals = [repeat['bulk_atoms_graphene'] for repeat in bulk_particles]
    bulk_graphene_avg_err = [np.average(bulk_graphene_vals), np.std(bulk_graphene_vals)]
    bulk_graphene_lost = bulk_atoms_dict["graphene"] - bulk_graphene_avg_err[0]

    penetration_depths = [particle.max_penetration for particle in particles]
    diamond_penetration_depths = [particle.diamond_max_penetration for particle in particles]
    avg_pen = np.average(penetration_depths)
    pen_err = np.std(penetration_depths)
    diamond_avg_pen = np.average(diamond_penetration_depths)
    diamond_pen_err = np.std(diamond_penetration_depths)

    finish_diamond = len([particle.finish for particle in particles if particle.finish == 'diamond'])
    finish_graphene = len([particle.finish for particle in particles if particle.finish == 'graphene'])
    finish_outside = len([particle.finish for particle in particles if particle.finish == 'outside'])

    pen_graphene = len([particle.pen for particle in particles if particle.pen == 'graphene' or particle.pen == 'graphene and diamond'])
    pen_diamond = len([particle.pen for particle in particles if particle.pen == 'graphene and diamond'])

    #output string for results text file
    output = "Bombardment Atom Energy: %seV"%round(energy,4)

    output += "\nSurface at: %s"%surface
    output += "\nNumber of graphene sheets: %s"%no_of_sheets
    output += "\nBombarding particle: %s"%atom_type
    output += "\nDiamond form: %s"%diamond_type
    output += "\nAverage Penetration: %s ± %s"%(round(avg_pen,6), round(pen_err,6))
    output += "\nAverage Diamond Penetration: %s ± %s"%(round(diamond_avg_pen,6), round(diamond_pen_err,6))
    output += "\nMaximum Penetration: %s, in repeat %s."%(round(max(penetration_depths),6),penetration_depths.index(max(penetration_depths)) + 1)
    output += "\nMaximum Diamond Penetration: %s, in repeat %s."%(round(max(diamond_penetration_depths),6),diamond_penetration_depths.index(max(diamond_penetration_depths)) + 1)
    output += '\nOutside material as final location: %s%% (%s particles).'%(round(finish_outside*100/len(particles),1), finish_outside)
    output += '\nGraphene as final location: %s%% (%s particles).'%(round(finish_graphene*100/len(particles),1), finish_graphene)
    output += '\nDiamond as final location: %s%% (%s particles).'%(round(finish_diamond*100/len(particles),1), finish_diamond)
    output += '\nPenetrated graphene: %s%% (%s particles).'%(round(pen_graphene*100/len(particles),1), pen_graphene)
    output += '\nPenetrated graphene and diamond: %s%% (%s particles).'%(round(pen_diamond*100/len(particles),1), pen_diamond)
    output += "\nAverage bulk diamond atoms displaced: {:3f} ± {:3f} out of {:}, {:3f}%%.".format(bulk_diamond_lost, bulk_diamond_avg_err[1],
                                                                                                 bulk_atoms_dict["diamond"], 
                                                                                                 (bulk_diamond_lost*100/bulk_atoms_dict["diamond"]))
    output += "\nAverage bulk graphene atoms displaced: {:3f} ± {:3f} out of {:}, {:3f}%%.".format(bulk_graphene_lost, bulk_graphene_avg_err[1],
                                                                                                bulk_atoms_dict["graphene"], 
                                                                                                (bulk_graphene_lost*100/bulk_atoms_dict["graphene"]))
                                                                                            

    output += '\n%s repeats found.\n\n'%len(particles)

    for particle in particles:
        output += particle.string

    with open("%s/penetration_results.txt"%path, 'w') as fp: 
        fp.write(output) 

    #results dict is used to move data into results_plotter.py 
    results = dict(energy = energy, atom_type = atom_type, diamond_type = diamond_type, surface = surface, 
                no_of_sheets = no_of_sheets, avg_pen = [avg_pen,pen_err], 
                diamond_avg_pen = [diamond_avg_pen,diamond_pen_err], max_pen = max(penetration_depths), 
                diamond_max_pen = max(diamond_penetration_depths), locations = [finish_diamond, finish_graphene, finish_outside], 
                pens = [pen_diamond, pen_graphene], repeats = len(particles), 
                diamond_atoms_lost = [value*100/bulk_atoms_dict["diamond"] for value in bulk_diamond_avg_err],
                graphene_atoms_lost = [value*100/bulk_atoms_dict["graphene"] for value in bulk_graphene_avg_err])

    return results

if __name__ == "__main__":
    path = '/Users/jamespittard/documents/culham/lammpsdepo/lammps/examples/bombardment/309.5_all_0'
    main(path, 1)

