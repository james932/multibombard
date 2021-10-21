



import sys
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd

dir_path = os.path.dirname(os.path.realpath(__file__))

files = sys.argv[1:]



pre_bomb_run_time = 3000
between_bomb_run_time = 3000
dump = 1000
last_t_step = 1506000
last_t_step = 4521000

def time_to_atom_conv(timestep, pre_bomb_run_time, between_bomb_run_time):

    bomb_time = timestep - pre_bomb_run_time

    if bomb_time > 0:
        bombs = int(bomb_time/between_bomb_run_time) + 1

    else:
        bombs = 0

    return bombs

def sort_name(path):
    timestep = path.split('/')[-1]
    timestep = float(timestep[:-4])
    return timestep

for file_count, file_ in enumerate(files):

    xyz_files = glob.glob(f"{dir_path}/hpc_results/{file_}*.xyz")
    xyz_files.sort(key = sort_name)

last_t_step_path = xyz_files[-1].split("/")
last_t_step = int(last_t_step_path[-1][:-4])



data = np.zeros((len(files),((last_t_step//between_bomb_run_time) +1),5))

for file_count, file_ in enumerate(files):

    xyz_files = glob.glob(f"{dir_path}/hpc_results/{file_}*.xyz")

    xyz_files.sort(key = sort_name)

    array_index = 0
    for index, path in enumerate(xyz_files):

        xyz = open(path, 'r')
        xyz = xyz.read()

        lines = xyz.split("\n")
        lines.remove('')

        timestep = path.split('/')[-1]
        timestep = float(timestep[:-4])


        timestep = float(lines[1])

        if timestep%between_bomb_run_time == 0:
            print(f"TIME STEP: {timestep}")
            atoms = float(lines[3])

            d_counter = 0
            t_counter = 0
            bulk_counter = 0

            for line in lines:
                line = line.split(' ')

                try:
                    if line[1] == '1':
                        bulk_counter += 1
                    if line[1] == '2':
                        d_counter += 1
                    if line[1] == '3':
                        t_counter += 1
                except IndexError:
                    pass

            
            bombs = time_to_atom_conv(timestep, pre_bomb_run_time, between_bomb_run_time)

            data[file_count][array_index][:][:][:] = [timestep, bombs, d_counter, t_counter, bulk_counter]
            array_index += 1
                
    def summative(list_, index):

        list_[index] = sum(list_[:index])
        return list_


for file_count,sim in enumerate(data):
    df = pd.DataFrame(sim)
    filepath = f"{dir_path}/hpc_results/data_array_{sys.argv[(1 + file_count)][:-1]}.csv"
    print(filepath)
    #with pd.ExcelWriter(filepath) as writer:
        #df.to_excel(writer, sheet_name='Sheet1')
    df.to_csv(filepath, index = False)


for sim in data:
    xs = [dat[1] for dat in sim if dat[1] != 0]
    bulk_raw = [dat[4] for dat in sim if dat[1] != 0]

    plt.plot(xs,bulk_raw)

plt.ylabel("Sample atoms")
plt.xlabel("Number of Bombarded Atoms")
plt.savefig(f"{dir_path}/hpc_results/{file_}d_perc.png")
plt.show()



for sim in data:
    xs = [dat[1]*4 for dat in sim if dat[1] != 0]
    d_perc = [dat[2]*100/dat[1] for dat in sim if dat[1] != 0]

    plt.plot(xs,d_perc)

plt.ylabel("Percentage of Bombarding Atoms Retained")
plt.xlabel("Number of Bombarded Atoms (adjusted)")
plt.ylim(-5,105)
plt.legend(["Diamond", "Diamond w/ Graphene"])
plt.savefig(f"{dir_path}/hpc_results/{file_}d_perc.png")
plt.show()

for sim in data:
    xs = [dat[1] for dat in sim if dat[1] != 0]
    d_raw= [dat[2] for dat in sim if dat[1] != 0]

    plt.plot(xs,d_raw)

plt.ylabel("Deuterium Atom Counts")
plt.xlabel("Number of Bombarded Atoms")
plt.legend(["Diamond", "Diamond w/ Graphene"])
plt.savefig(f"{dir_path}/hpc_results/{file_}d_raw.png")
plt.show()


'''         
d_xs = [dat[1] for dat in data if dat[1] != 0]
#d_ys = [dat[2]*100/data[0][2] for dat in data if dat[1] != 0] #trit
#d_ys = [dat[2]*100/(dat[1]+ 2000) for dat in data if dat[1] != 0] #dcont
d_ys = [dat[2]*100/dat[1] for dat in data if dat[1] != 0] #d

t_xs = [dat[1] for dat in data if dat[1] != 0]
t_ys = [dat[3]*100/dat[1] for dat in data if dat[1] != 0]

raw_ys = [dat[2] for dat in data if dat[1] != 0]
raw_ts = [dat[3] for dat in data if dat[1] != 0]

bulk = [dat[4] for dat in data]
bulk_xs = [dat[1] for dat in data]



plt.plot(bulk_xs[:-10], bulk[:-10])
plt.title(f"{file_}")
plt.ylabel("Bulk Atoms")
plt.xlabel("Time / ps")
plt.savefig(f"{dir_path}/hpc_results/{file_}bulk.png")
plt.show()

plt.plot(d_xs,d_ys)
plt.title(f"{file_}")
plt.ylabel("D Percentage Retained")
plt.xlabel("Number of bombarded atoms")
plt.ylim(-5,105)
plt.savefig(f"{dir_path}/hpc_results/{file_}d_perc.png")
plt.show()

plt.plot(t_xs,t_ys)
plt.title(f"{file_}")
plt.ylabel("T Percentage Retained")
plt.xlabel("Number of bombarded atoms")
plt.ylim(-5,105)
plt.savefig(f"{dir_path}/hpc_results/{file_}t_perc.png")
plt.show()

plt.plot(d_xs, raw_ys)
plt.title(f"{file_}")
plt.ylabel("D Number Retained")
plt.xlabel("Number of bombarded atoms")
plt.savefig(f"{dir_path}/hpc_results/{file_}d_num.png")
plt.show()

plt.plot(d_xs, raw_ts)
plt.title(f"{file_}")
plt.ylabel("T Number Retained")
plt.xlabel("Number of bombarded atoms")
plt.savefig(f"{dir_path}/hpc_results/{file_}t_num.png")
plt.show()
'''