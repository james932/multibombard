


import numpy as np
from collections import namedtuple


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

            

class Frame():

    def __init__(self, time, atoms):
        self.time  = time
        self.atoms = atoms

        #sortAtoms(self)


    def __getitem__(self, index):
        out = None
        for atom in self.atoms:
            if atom.index == float(index):
                out = atom
        return out


    def sortAtoms(self):

        atoms = self.atoms
        sorted_atoms = []

        for i in range(len(atoms)):
            for atom in atoms:
                if i == atom.index:
                    sorted_atoms.append(atom)
                    atoms.remove(atom)
                    

        setattr(self, 'atoms', sorted_atoms)
        return sorted_atoms

    def present(self, jmol = False):
        frame = ''
        frame += '%s'%len(self.atoms)
        frame += '\nAtoms. Timestep. %s'%self.time
        
        if jmol == True:
            for atom in self.atoms:
                frame += f"\n{atom.atom_type} {atom.pos[0]} {atom.pos[1]} {atom.pos[2]}"

        else:
            for atom in self.atoms:
                frame += f"\n{atom.index} {atom.atom_type} {atom.pos[0]} {atom.pos[1]} {atom.pos[2]}"

        frame += '\n\n'

        return frame



class Results():

    def __init__(self, first_frame, last_frame, all_frames = None, no_of_sheets = None, 
                atom_type = None, energy = None, diamond_type = None, replicate = None, 
                no_bombarding_atoms = None):

        details = namedtuple('details', ['no_of_sheets', 'atom_type', 'energy', 
                                        'diamond_type', 'replicate', 'no_bombarding_atoms'])
        
        first_frame.sortAtoms()
        last_frame.sortAtoms()

        self.first_frame = first_frame
        self.last_frame = last_frame
        self.all_frames = all_frames
        self.details = details(int(no_of_sheets), atom_type, energy, diamond_type, 
                                replicate, no_bombarding_atoms)

        self.surface = (None, None)

    def getRetainedParticles(self):
        retained_particles = [atom for atom in self.last_frame.atoms if atom.atom_type == 2]
        setattr(self, 'retained_particles', retained_particles)
        return len(retained_particles)

    def findRegions(self):

        first_frame_atoms = self.first_frame.atoms
        diamond_atoms = [atom for atom in first_frame_atoms if atom.pos[2] >= 0]
        diamond_surface_atoms = [atom for atom in first_frame_atoms if atom.pos[2] == 0]
        unsorted_graphene_atoms = [atom for atom in first_frame_atoms if atom.pos[2] < 0]
        
        graphene_atoms = [[] for i in range(self.details.no_of_sheets)]
        print(f"{len(unsorted_graphene_atoms)} graphene atoms in total.")
        for i in range(self.details.no_of_sheets):
            graphene_atoms[i] = [atom for atom in unsorted_graphene_atoms if atom.pos[2] == round((-1.6775 - i*3.35),6)]
            print(f"Graphene sheet {i} has {len(graphene_atoms[i])} atoms and a height of {-1.6775 - i*3.35}.")


        setattr(self, 'diamond_atoms', diamond_atoms)
        setattr(self, 'graphene_atoms', graphene_atoms)

        last_frame_atoms = self.last_frame.atoms 
        
        diamond_surface_indexes = [atom.index for atom in diamond_surface_atoms]
        
        diamond_surface_heights = [atom.pos[2] for atom in last_frame_atoms if atom.index in diamond_surface_indexes]
        diamond_surface = (np.average(diamond_surface_heights), np.std(diamond_surface_heights))

        layer_surfaces = []
        for layer_atoms in graphene_atoms:
            indexes = [atom.index for atom in layer_atoms]

            layer_heights = [atom.pos[2] for atom in last_frame_atoms if atom.index in indexes]                 
            layer_surfaces.append((np.average(layer_heights), np.std(layer_heights)))
            

        setattr(self, "diamond_surface", diamond_surface)
        setattr(self, "graphene_surfaces", layer_surfaces)
        

        print(f"diamond surface: {diamond_surface[0]} ± {diamond_surface[1]}")
        for layer in layer_surfaces:
            print(f"graphene surface: {layer[0]} ± {layer[1]}")

        
        

    def getAtomMovement(self):

        displacement_mag = []
        
        for i in range(len(self.first_frame.atoms)):

            try:
            
                if self.last_frame.atoms[i].atom_type == 1:
                    r = self.last_frame.atoms[i].pos - self.first_frame.atoms[i].pos

                    r_abs = np.sqrt(r.dot(r))
                     
                    if r_abs < 10:
                        displacement_mag.append(r_abs)

            except IndexError:
                r_abs = None


        setattr(self, 'displacement_mag', displacement_mag)
        return max(displacement_mag)

    def analyseLastPositions(self):

        #retained_atoms = [atom for atom in self.last_frame.atoms if atom.atom_type == 2]
        surface_heights = [self.diamond_surface[0]] + [layer[0] for layer in self.graphene_surfaces]
        surface_heights.sort(reverse = True)

        diamond_atoms = []
        graphene_atoms = [[] for i in range(self.details.no_of_sheets + 1)]

        for atom in self.retained_particles:
    
            surface_heights.append(atom.pos[2])
            surface_heights.sort(reverse = True)
            index = surface_heights.index(atom.pos[2])
            surface_heights.remove(atom.pos[2])
            if index == 0:
                diamond_atoms.append(atom)

            if index > 0:
                graphene_atoms[index-1].append(atom)

                
        diamond_pens = [atom.pos[2] - self.diamond_surface[0] for atom in diamond_atoms]
        avg_diamond_pen = (np.average(diamond_pens), np.std(diamond_pens))

        surface_pens = [atom.pos[2] - min(surface_heights) for atom in self.retained_particles]
        avg_surface_pen = (np.average(surface_pens), np.std(surface_pens))

        number_penetrating_graphene = [len(layer) for layer in graphene_atoms]
        number_penetrating_diamond = len(diamond_atoms)

        last_frame_info = namedtuple('last_frame_info', ['retained_particles','diamond_atoms', 'graphene_atoms',
                                                        'avg_diamond_pen', 'avg_surface_pen'])
                                            
        info = last_frame_info(len(self.retained_particles), len(diamond_atoms), number_penetrating_graphene, 
                                avg_diamond_pen, avg_surface_pen)

        setattr(self, 'last_frame_info', info)

    def getDepthProfile(self):

        if self.diamond_surface[0] != None:
            surface = self.diamond_surface[0]
        else:
            surface = 0

        diamond_depths = [atom.pos[2] - surface for atom in self.retained_particles]
        pen_diamond_depths = [depth for depth in diamond_depths if depth >0]
        setattr(self, 'diamond_depth_profile', pen_diamond_depths)

        if len(self.graphene_surfaces) != 0:
            surfaces = [surface[0] for surface in self.graphene_surfaces]
            surface = min(surfaces)

            depths = [atom.pos[2] - surface for atom in self.retained_particles]

        else:
            depths = diamond_depths

        setattr(self, 'surface_depth_profile', depths)
          

    def present(self):

        counter = 0
        while True:
            #try:
            ##    if counter > 10:
             #       break

            out = '\n-------RESULTS--------\n'
            out += "\nBombarding Atom Type: %s"%self.details.atom_type
            out += "\nBombarding Atom Energy: %seV"%self.details.energy
            out += "\nDiamond Type: %s"%self.details.diamond_type
            out += "\nNumber of Graphene Layers %s"%self.details.no_of_sheets
            out += "\nMaximum Bulk Atom Displacement: %s"%max(self.displacement_mag)
            out += f"\nDiamond Surface height: {self.diamond_surface[0]} ± {self.diamond_surface[1]}"
            for i in range(self.details.no_of_sheets):
                out += f"\nGraphene Layer {i+1} height: {self.graphene_surfaces[i][0]} ± {self.graphene_surfaces[i][1]}"
            out += "\nNumber of Bulk Atoms: %s"%(len(self.first_frame.atoms) - 1) 
            out += f"\nNumber of Bombardment Atoms: {self.details.no_bombarding_atoms}"
            out += f"\nNumber of retained bombardment atoms: {self.last_frame_info.retained_particles}"
            out += f"\nNumber of bombardment atoms in diamond: {self.last_frame_info.diamond_atoms}"
            out += f"\nNumber of bombardment atoms in graphene: {self.last_frame_info.graphene_atoms}"
            out += f"\nNote above refers to regions between diamond and layer 1, layer 1 and layer 2 etc."
            out += f"\nAverage diamond penetration: {self.last_frame_info.avg_diamond_pen[0]} ± {self.last_frame_info.avg_diamond_pen[1]}"
            out += f"\nAverage surface penetration: {self.last_frame_info.avg_surface_pen[0]} ± {self.last_frame_info.avg_surface_pen[1]}"
            out += f"\n\nSurface Depth Profile:"
            out += f"\n{self.surface_depth_profile}"
            out += f"\n\nDiamond Depth Profile:"
            out += f"\n{self.diamond_depth_profile}"
            out += '\n'
            break
            '''
            except AttributeError:
               
                self.getRetainedParticles()
                self.getAtomMovement()
                self.findRegions()
                self.analyseLastPositions()
                self.getDepthProfile()
                counter += 1
            '''
        return out


    


def main(all_xyz_file_path, no_of_sheets, atom_type, energy, diamond_type, replicate, no_bombarding_atoms):
        # , no_of_sheets, atom_type, energy, diamond_type, surface = 0, bulk_atoms_dict = None):
    
    atom = namedtuple('atom', ['index', 'atom_type', 'pos', 'time'])

    frames = []

    all_xyz = open(all_xyz_file_path, 'r')
    all_xyz = all_xyz.read()

    all_xyz = all_xyz.split("ITEM: TIMESTEP")
    all_xyz.remove('')

    jmol_str = ''
    
    #for time_step in all_xyz:
    for i in range(0,len(all_xyz), 100):

        time_step = all_xyz[i]

        time_step = time_step.split('\n')
        time_step.remove('')
 
        time = float(time_step[0])
        atoms = []

        jmol_atom_str = ''

        for i in time_step:
            line = i.split()

            if len(line) == 5:
                index = float(line[0])
                atom_type_no = float(line[1])
                atom_pos = np.array([float(line[2]), float(line[3]), float(line[4])])

                atoms.append(atom(index, atom_type_no, atom_pos, time))
                jmol_atom_str += "\n%s %s %s %s"%(line[1], line[2], line[3], line[4])

        frame = ''
        frame += '%s'%len(atoms)
        frame += '\nAtoms. Timestep. %s'%time
        frame += jmol_atom_str
        frame += '\n\n'
        jmol_str += frame

        frames.append(Frame(time, atoms))


    with open("%s/jmol_all.xyz"%(all_xyz_file_path[:-7]), 'w') as fp: #rewriting edited input file
        fp.write(jmol_str)

    with open("%s/initial.xyz"%(all_xyz_file_path[:-7]), 'w') as fp: #rewriting edited input file
        fp.write(frames[0].present(jmol = True))

    with open("%s/final.xyz"%(all_xyz_file_path[:-7]), 'w') as fp: #rewriting edited input file
        fp.write(frame)

    with open("%s/final_test.xyz"%(all_xyz_file_path[:-7]), 'w') as fp: #rewriting edited input file
        fp.write(frames[-1].present())

    

    results = Results(first_frame=frames[0], last_frame=frames[-1], all_frames = None, 
                        no_of_sheets = no_of_sheets, atom_type =  atom_type, energy = energy,
                        diamond_type = diamond_type, replicate = replicate, 
                        no_bombarding_atoms = no_bombarding_atoms)
    
    results.getAtomMovement()
    results.getRetainedParticles()
    results.findRegions()
    results.analyseLastPositions()
    results.getDepthProfile()

    return results


