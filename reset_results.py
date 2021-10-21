





import os


dir_path = os.path.dirname(os.path.realpath(__file__))
results_path = dir_path + "/results/"

os.system(f"rm -r {dir_path}/results/d_0g_80eV_11/")
os.system(f"mkdir {dir_path}/results/d_0g_80eV_11/")

os.system(f"cp -r {dir_path}/results/sample_res/* {dir_path}/results/d_0g_80eV_11/")