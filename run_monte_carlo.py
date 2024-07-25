from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import math
from footprint_filter import filter_by_footprint
from voidiness import voidy_analysis
import os # Validating user input

# Defaults
DEF_CAT_FN = 'exported_dataFrames/4lac_w_voidiness.xlsx'
DEF_VOID_FN = 'exported_dataFrames/voids.xlsx'
DEF_FP_FN = 'exported_dataFrames/footprint_points_void_centers.xlsx'
MEM_LIM = 10 # Memory limit in megabytes
def get_usr_in():
    print("Default celestial object catalog: " + DEF_CAT_FN)
    print("Default Void:                     " + DEF_VOID_FN)
    print("Default footprint:                " + DEF_FP_FN)

    ans = input("Would you like to run a different CELESTIAL OBJECT catalog? (Only excel files accepted) [y/(n)] ")

    assert (ans == "y" or ans == "n" )or ans=='', f"Must type 'y' or 'n'. Recieved: {ans}"
    if ans == 'y':
        catalog_fn = input("Type in path to file name, relative or absolute:")
        assert os.path.isfile(catalog_fn), f"File does not exist: {catalog_fn}"
    elif ans =='n' or ans == '':
        catalog_fn = DEF_CAT_FN

    ans = input("Would you like to run a different VOID catalog? (Only excel files accepted) [y/(n)] ")

    assert (ans == "y" or ans == "n" )or ans=='', f"Must type 'y' or 'n'. Recieved: {ans}"
    if ans == 'y':
        void_fn= input("Type in path to file name, relative or absolute:")
        assert os.path.isfile(catalog_fn), f"File does not exist: {void_fn}"
    elif ans =='n' or ans == '':
        void_fn = DEF_VOID_FN

    ans = input("Would you like to run a different FOOTPRINT? (Only excel files accepted) [y/(n)] ")
    assert (ans == "y" or ans == "n" )or ans=='', f"Must type 'y' or 'n'. Recieved: {ans}"
    if ans == 'y':
        fp_fn= input("Type in path to file name, relative or absolute:")
        assert os.path.isfile(catalog_fn), f"File does not exist: {fp_fn}"
    elif ans =='n' or ans == '':
        fp_fn = DEF_FP_FN
    
    print("\nUsing:\n",
          "Celestial object catalog: " + catalog_fn + "\n",
          "Void catalog:             " + void_fn+ "\n",
          "Footprint:                " + fp_fn +  "\n")
    
    ans = input(f"Change memory limit? Current: {MEM_LIM} MB (y/[n])\n")
    assert (ans == "y" or ans == "n" ) or ans=='', f"Must type 'y' or 'n'. Recieved: {ans}"
    if ans == "y":
        mem_lim = float(input("Input new memory limit in MB: "))
    else:
        mem_lim = MEM_LIM

    return (catalog_fn, void_fn, fp_fn, mem_lim)

def rand_long_and_lat(n, seeded=False, seed = 567307250717):
    """
    Generates number of longitude and latitude coordinates in degrees for entire pandasDF
    """
    if seeded:
        seed = seed # Using seed while debugging
    else:
        seed = None

    rng = np.random.default_rng(seed) 
    theta = np.arccos(1 - 2 * rng.uniform(0, 1, n)) * (180/math.pi) # COLAT
    b = 90 - theta
    l = rng.uniform(0, 360, n)
    return (l, b)
# Randomize

def gen_filtered_ra_deg(n, footprint_fn):
    randRA, randDE = rand_long_and_lat(n)
    temp_cel = pd.DataFrame({'RAdeg': randRA, 'DEdeg': randDE})
    temp_cel = filter_by_footprint(temp_cel, fp_fn)
    return temp_cel

def monte_carlo(cat_fn, void_fn, fp_fn):
    # voids = pd.read_excel('exported_dataFrames/voids.xlsx')
    voids = pd.read_excel(void_fn)
    master_Carlo = pd.read_excel(cat_fn)

    run_n = 0
    mc_voidiness = np.array([])
    n_galxy = len(master_Carlo)

    while True:
        coords = gen_filtered_ra_deg(100000)
        coords_idx = coords.index.tolist()
        while len(coords_idx) > n_galxy:
            fresh_coords_idx = coords_idx[:n_galxy]
            coords_idx = coords_idx[n_galxy:]

            master_Carlo['RAdeg'] = coords.RAdeg[fresh_coords_idx].values
            master_Carlo['DEdeg'] = coords.DEdeg[fresh_coords_idx].values
            mc_voidiness = np.append(mc_voidiness, voidy_analysis(voids, master_Carlo).Voidiness.values) # This one takes about 6s per run
            run_n += 1
   

if __name__ == "__main__":
    print(get_usr_in())
    # cat_fn, void_fn, fp_fn, mem_lim = get_usr_in()
    # monte_carlo(cat_fn, void_fn, fp_fn)