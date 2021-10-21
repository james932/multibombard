

import os

def file_proc(file):
    """Short function to split text files into a list of lines"""
    opened_file = open(file, 'r')
    opened_file = opened_file.read()
    opened_file = opened_file.split("\n")
    return opened_file

def atom_locator(diamond_data_file):

    dir_path = os.path.dirname(os.path.realpath(__file__))
    og_data_file = file_proc("%s/LAMMPS_files/%s"%(dir_path,diamond_data_file))

    data = False
    atoms = []
    for line in og_data_file:
        
        if line != '':
            if data == True:

                coord_line = line.split()
                coord = coord_line[2:]
                coord = [float(i) for i in coord]

                atoms.append(coord)

        if line == "Atoms":
            data = True

    xs = [pos[0] for pos in atoms]
    ys = [pos[1] for pos in atoms]
    zs = [pos[2] for pos in atoms]

    diamond_dims = dict(x_lims = [min(xs), max(xs)], y_lims = [min(ys), max(ys)],
                        z_lims = [min(zs), max(zs)])

    return diamond_dims



class Sheet:
    """Sheet object to hold atoms contained within one sheet of graphene. """
    def __init__(self,number, atoms):

        self.number = number 
        self.atoms = atoms
        if number != None:
            self.height = atoms[0].z


class Atom:
    """Atom object used ofr each carbon atom in graphene. """
    def __init__(self, atom):

        number, mass_index, position = atom
        self.mass_index = str(mass_index)
        self.id = str(number)
        self.mass = 1
        self.x = float(position[0])
        self.y = float(position[1])
        self.z = float(position[2])
        self.position = str("%s %s %s " %(position[0], position[1],position[2]))

        self.string = "%s %s %s"%(self.id, self.mass_index, self.position)



def main(path = None, data_file_name = None, input_file_name = None, new_data_file_name = None, atom_mass = None):
    """Creates graphene data file of number of sheets specified in the input file."""

    if __name__ == "__main__": #Can be used to creat data file outside of the queue package
        path = "/Users/jamespittard/Documents/Culham/lammpsdepo/lammps/examples/diamond_graphene/"
        data_file_name = 'data.graphene'

    in_file = file_proc("%s/%s"%(path,input_file_name))

    if data_file_name != None:

        data_file = file_proc("%s/%s"%(path,data_file_name))

        atom_switch = False #switch used to identify where data needs to be extracted
        unit_atoms = []
        for row in data_file: #goes through rows in data files

            if row == 'Atoms':
                atom_switch = True
                atom_index = data_file.index(row)
            try:
                if row.split()[2] == "xlo": #box parameters
                    x_para = row.split()[0:2]
                if row.split()[2] == "ylo":
                    y_para = row.split()[0:2]
                if row.split()[2] == "zlo":
                    z_para = row.split()[0:2]
                
                if atom_switch == True: #creates atom objects
                    row = row.split()
                    row = [row[0], row[1], row[2:]]
                    atom = Atom(row)
                    unit_atoms.append(atom)

            
            except IndexError:
                pass

        #box lengths
        x_len = abs(float(x_para[0]) - float(x_para[1]))
        y_len = abs(float(y_para[0]) - float(y_para[1]))
        z_len = abs(float(z_para[0]) - float(z_para[1]))

    replicate = None

    for row in in_file: #reading input file to determine size and number of sheets
        
        row = row.split()
        try:
            if row[0] == "variable":
                no_bombarding_atoms = float(row[-1])

            if row[0] == 'replicate':
                replicate = row[1:]
            
            if row[0] == '#number_of_sheets':
                no_of_sheets = float(row[1])

            if row[0] == "#diamond_type":
                diamond_type = row[1]
            
            if row[0] == "read_data":
                diamond_path = row[1].split("/")
                if diamond_path[-1][:12] == "data.diamond":
                    diamond_data_file = diamond_path[-1]
                    #this is used if replicate is not being used
                    #but rather a larger data file with all atom postions

            
            if row[0] == "#atom_type":
                atom_type = row[1]
                print(f"gmak atom_type: {atom_type}")
        except IndexError:
            pass

    if replicate == None: #one larger data file used for diamond instead
        replicate = [8, 8, 10]
        diamond_dims = atom_locator(diamond_data_file)
        xlo, xhi = diamond_dims["x_lims"]
        ylo, yhi = diamond_dims["y_lims"]
        z_shift = diamond_dims["z_lims"][0]
    else:
        replicate = [int(i) for i in replicate]
        diamond_size = 3.567
        xhi = replicate[0]*diamond_size
        yhi = replicate[1]*diamond_size
        z_shift = 0

    if data_file_name != None:
        count = 1
        new_atoms = []
        sheet_atoms = []
        diamond_size = 3.567
        graphene_thickness = 3.35

        if no_of_sheets == 0:
            new_atoms = []
            origin_sheet = Sheet(None, [])

        else: #creates first sheet spanning surface
            for x_shift in range(0,int(replicate[0])):
                for y_shift in range(0, int(replicate[1])):

                    for atom in unit_atoms: #adds on unit cell shifts to original data file
                        shift = [x_len*x_shift, y_len*y_shift, z_shift]
                        new_position = [round(atom.x + shift[0],6), round(atom.y + shift[1],6), round(-atom.z + shift[2],6)] 
                    
                        mass_index = atom.mass_index

                        if new_position[0] < xhi and new_position[1] < yhi:
                            new_atom = Atom([count,mass_index,new_position])
                            sheet_atoms.append(new_atom)
                            count +=1

            origin_sheet = Sheet(0,sheet_atoms)

            sheet_no = 1
            sheets = [origin_sheet]

            while sheet_no < no_of_sheets: #creates further sheets from origin sheet

                sign = 1 #used to add offset to consecutive layers
                sheet_atoms = []

                if sheet_no%2 == 0:
                    sign = 0

                for atom in origin_sheet.atoms:
                    new_position = [round(atom.x + sign*1.42,6), round(atom.y,6), round(atom.z - sheet_no*3.35,6)] 
                    
                    mass_index = atom.mass_index

                    new_atom = Atom([count,mass_index,new_position])
                    sheet_atoms.append(new_atom)
                    count +=1

                sheets.append(Sheet(sheet_no, sheet_atoms))
                sheet_no += 1

            for sheet in sheets:
                new_atoms += sheet.atoms


        number_of_atoms = len(new_atoms)
        data_file_start = data_file[:atom_index+1] #selecting data file up to the atom data 
        atoms_string = ''

        for item in data_file_start: #reads and edits start of data file
            try:
                if item.split()[1] == 'atoms': 
                    item = '%s atoms'%number_of_atoms
            
                if item.split()[0] == '2' and len(item.split()) == 2:
                    item = "2 %s"%atom_mass

                if item.split()[2] == "xlo":
                    item = "0 %a xlo xhi"%(xhi)        
                if item.split()[2] == "ylo":
                    item = "0 %s ylo yhi"%(int(replicate[1])*diamond_size)
                if item.split()[2] == "zlo":
                    item = "%s %s zlo zhi"%((-no_of_sheets*graphene_thickness*2),(int(replicate[2])*diamond_size*1.5))

            except IndexError:
                pass

            atoms_string +=  item + "\n" 


        for atom in new_atoms: # adding atom coordinates
            atoms_string += '\n' + atom.string

        with open ('%s/%s'%(path, new_data_file_name), 'w') as fp: #writing new data file
            print('%s/%s'%(path, new_data_file_name))
            fp.write(atoms_string) 

    #used later, gives the initial number of atoms in each region.
    if data_file_name == None:

        data_file = file_proc("%s/data.graphite_sheet"%path)

        atoms_line = data_file[2].split()
  
        number_of_atoms = int(atoms_line[0])


    no_of_bulk_atoms = dict(diamond = 8*int(replicate[0])*int(replicate[1])*int(replicate[2]), graphene = number_of_atoms)

    if data_file_name == None:
        return no_of_bulk_atoms, no_of_sheets, atom_type, diamond_type, replicate, no_bombarding_atoms

    else:
        return no_of_bulk_atoms


