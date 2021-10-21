


import os
import sys
import numpy



file_ = None


dir_path = os.path.dirname(os.path.realpath(__file__))

file_ = sys.argv[1]


file_ = open(f"{dir_path}/hpc_results/{file_}", 'r')
file_ = file_.read()

file_ = file_.split("\n")
file_.remove('')

atom_str = ''

for index, line in enumerate(file_):
    line = line.split(' ')

    if index == 3:
        atoms = int(line[0])

    if index == 5:
        xlo, xhi = line

    if index == 6:
        ylo,yhi = line

    if index == 7:
        zlo,zhi = line

    if index > 8:
        atom_str += f"\n{index -8} {line[1]} {line[2]} {line[3]} {line[4]}"

data_file = 'LAMMPS data file from restart file: timestep = 1, procs = 1'
data_file += f'\n\n{atoms} atoms'
data_file += "\n\n3 atom types"
data_file += f"\n\n{xlo} {xhi} xlo xhi"
data_file += f"\n{ylo} {yhi} ylo yhi"
data_file += f"\n{zlo} {zhi} zlo zhi"
data_file += "\n\nMasses"
data_file += "\n\n1 12.01"
data_file += "\n2 2.0014"
data_file += "\n3 3.0160" 
data_file += "\n\n Atoms"
data_file += "\n" + atom_str



with open (('%s/LAMMPS_files/data.new_data'%dir_path), 'w') as fp: #writing new data file
    fp.write(data_file) 
