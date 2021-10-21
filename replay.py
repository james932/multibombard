




import os
import numpy as np
import sys
import matplotlib.pyplot as plt


#to be run from the results dir
dir_path = os.path.dirname(os.path.realpath(__file__))

pre_bomb_run_time = 3000
between_bomb_run_time = 3000
dump = 1000

file_ = sys.argv[1]


all_file = open(f"{dir_path}/hpc_results/{file_}all.xyz", 'r')
all_file = all_file.read()

#target_number = sys.argv[2]


timesteps = all_file.split("ITEM: TIMESTEP")
timesteps.remove('')

data = np.zeros((len(timesteps),4))


def time_to_atom_conv(timestep, pre_bomb_run_time, between_bomb_run_time):

    bomb_time = timestep - pre_bomb_run_time

    if bomb_time > 0:
        bombs = int(bomb_time/between_bomb_run_time) + 1

    else:
        bombs = 0

    return bombs

for index, frame in enumerate(timesteps):
    lines = frame.split("\n")
    lines.remove('')
    print(f"TIME STEP: {lines[0]}")
    timestep = float(lines[0])
    atoms = float(lines[2])

    counter = 0
    for line in lines:
        line = line.split(' ')

        try:
            if line[1] == '2':
                counter += 1
        except IndexError:
            pass

    bulk_atoms = atoms - counter
    bombs = time_to_atom_conv(timestep, pre_bomb_run_time, between_bomb_run_time)

    data[index][:][:][:] = [timestep, bombs, counter, bulk_atoms]


print(f"Initially: {data[0]} ")
print(f"MIddlely: {data[int(data.shape[0]/2)]} ")
print(f"finally: {data[-1]}")

xs = [dat[1] for dat in data if dat[1] != 0]
ys = [dat[2]*100/dat[1] for dat in data if dat[1] != 0]
raw_ys = [dat[2] for dat in data if dat[1] != 0]
bulk = [dat[3] for dat in data]
bulk_xs = [dat[0] for dat in data]

plt.plot(bulk_xs, bulk)
plt.ylabel("Bulk Atoms")
plt.xlabel("Time / ps")
plt.savefig(f"{file_}bulk.png")
plt.show()

plt.plot(xs,ys)
plt.ylabel("Percentage Retained")
plt.xlabel("Number of bombarded atoms")
plt.savefig(f"{file_}perc.png")
plt.show()

plt.plot(xs, raw_ys)
plt.ylabel("Number retained Retained")
plt.xlabel("Number of bombarded atoms")
plt.savefig(f"{file_}num.png")
plt.show()

'''



final = open(f"{dir_path}/hpc_results/{file_}final.xyz")
final = final.read()

final = final.split("\n")
count = 0
for line in final:

    line = line.split(' ')

    if line[0] == '2':
        count += 1

print(f"File {file_} retained {count} atoms")
'''
   