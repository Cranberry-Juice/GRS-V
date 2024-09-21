from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import math
from footprint_filter import filter_by_footprint
from voidiness import voidy_analysis
import os # Validating user input
from sys import getsizeof
from multiprocessing import Pool # Used for multicore processing
import pickle
import re
import time # time the code
# Defaults
DEF_CAT_FN = 'exported_dataFrames/qsos_w_voidiness.xlsx'
DEF_VOID_FN = 'exported_dataFrames/voids.xlsx'
DEF_FP_FN = 'exported_dataFrames/footprint_points_void_centers.xlsx'
MEM_LIM = 200 # Memory limit in megabytes
DEF_N_SAMP = 4 # Default number of samplings. 
DEF_CORE = 2 # Default number of cores to use
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

    ans = input(f"Change number of samplings? ({DEF_N_SAMP}) (y/[n]): ")
    assert (ans == "y" or ans == "n" ) or ans=='', f"Must type 'y' or 'n'. Recieved: {ans}"
    if ans == "y":
        n_samples = int(input("Input new number of samples must be EVEN number (minimum 4): "))
        while n_samples%2 != 0 and n_samples<4:
            n_samples = int(input(f"Input new number of samples must be EVEN number. Received: {n_samples}\n"))
    else:
        n_samples = DEF_N_SAMP
    # Prompt for number of CPU cores to use
    cores = input(f"Change number of CPU cores to use? (default: {DEF_CORE}) (y/[n]): ")
    assert (cores == "y" or cores == "n" or cores == ''), f"Must type 'y' or 'n'. Received: {cores}"

    if cores == "y":
        num_cores = int(input("Input number of CPU cores to use (minimum 1): "))
        while num_cores < 1:
            num_cores = int(input(f"Input number of CPU cores to use (minimum 1). Received: {num_cores}\n"))
    else:
        num_cores = DEF_CORE  # default value if not changing

    return (catalog_fn, void_fn, fp_fn, mem_lim, n_samples, num_cores)

def rand_long_and_lat(n, seed):
    """
    Generates number of longitude and latitude coordinates in degrees for entire pandasDF
    """
    # if seeded:
    #     seed = seed # Using seed while debugging
    # else:
    #     seed = None

    rng = np.random.default_rng(seed) 
    theta = np.arccos(1 - 2 * rng.uniform(0, 1, n)) * (180/math.pi) # COLAT
    b = 90 - theta
    l = rng.uniform(0, 360, n)
    return (l, b)
# Randomize

def gen_filtered_ra_deg(n, footprint_fn, seed):
    randRA, randDE = rand_long_and_lat(n, seed = seed)
    temp_cel = pd.DataFrame({'RAdeg': randRA, 'DEdeg': randDE})
    temp_cel = filter_by_footprint(temp_cel, footprint_fn)
    return temp_cel

def monte_carlo(n_samples, cat_fn, void_fn, fp_fn, mem_lim, seed):
    # voids = pd.read_excel('exported_dataFrames/voids.xlsx')
    voids = pd.read_excel(void_fn)
    master_Carlo = pd.read_excel(cat_fn)

    run_n = 1
    # mc_voidiness = np.array([])
    mc_voidiness = pd.DataFrame()
    n_galxy = len(master_Carlo)

    sim_size = getsizeof(mc_voidiness)/1e6 # size of simulated data in MB

    too_big = sim_size > mem_lim

    n_rand_samples = int(n_galxy * 1.5)
    rng = np.random.default_rng(seed)
    microseeds = [ int(rng.uniform(0, seed)) for i in range(2000)]
    while ~too_big and run_n < n_samples:
        # n_rand_samples = 100000

        # if n_galxy > 15000:
        #     # Mostly for the sdss catalog since it is already 100,000 big
        #     n_rand_samples*=10

        coords = gen_filtered_ra_deg(n_rand_samples, footprint_fn=fp_fn, seed=microseeds.pop())
        # print(f"OG Coords shape {coords.shape}")
        while len(coords) < n_galxy * 1.5:
            # print(f"n_rand_samples too small, generating more points. Now: {n_rand_samples}")
            coords = gen_filtered_ra_deg(n_rand_samples, footprint_fn=fp_fn, seed=microseeds.pop())
            # print(f"New Coords shape {coords.shape}")
            n_rand_samples += n_rand_samples # update the typical necessary samples to not run into this problem. 
        # print(repr(coords.shape))
        coords_idx = coords.index.tolist()
        # print(len(coords))
        # print(n_galxy)
        t0 = time.time()
        while len(coords_idx) > n_galxy:
            fresh_coords_idx = coords_idx[:n_galxy]
            coords_idx = coords_idx[n_galxy:]
            # print(len(fresh_coords_idx))

            master_Carlo['RAdeg'] = coords.RAdeg[fresh_coords_idx].values
            master_Carlo['DEdeg'] = coords.DEdeg[fresh_coords_idx].values
            new_voidy_values =  voidy_analysis(voids, master_Carlo)
            temp = new_voidy_values[['z', 'Voidiness']].copy()
            # mc_voidiness = np.append(mc_voidiness, voidy_analysis(voids, master_Carlo).Voidiness.values) # This one takes about 6s per 
            mc_voidiness = pd.concat((mc_voidiness, temp))
            t1 = time.time()
            print(f"Appended sample number {run_n} in {t1-t0:.2f} seconds!")
            t0 = time.time()
            if run_n >= n_samples or getsizeof(mc_voidiness)/1e6 > mem_lim:
                break
            run_n += 1
        too_big = getsizeof(mc_voidiness)/1e6 > mem_lim

        if run_n %10 == 0:
            # Print status every 10 runs
            print("Ran out of samplings")
            print(run_n)

    foo = 1 if too_big else 0
    bar = 2 if run_n >= n_samples else 0
    exitcode = foo + bar

    
    return (mc_voidiness, exitcode)

def save_dat(cat_fn, data):
    # Regular expression to match the part between '/' and '.xlsx'
    pattern = r'\/(.*?)\.xlsx$'
    name = re.search(pattern, cat_fn)
    name = name.group(1)
    fn = "stats/simulated_data/" + name + "_simulated_data.pkl"

    if os.path.isfile(fn):
        with open(fn, 'rb') as f:
            old_list = pickle.load(f)
        assert isinstance(old_list, pd.DataFrame), "Did not load Pandas Dataframe"
        print(f"Old length: {len(old_list)}")
        print(f"Adding {len(data)} data points to saved data")
        # data = np.append(old_list, data)
        data = pd.concat((old_list, data)).copy()

    with open(fn, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"File saved to: {fn}")
    print(f"{len(data)} total saved data points")

if __name__ == "__main__":
    cat_fn, void_fn, fp_fn, mem_lim, n_samples, n_core = get_usr_in()

    # Generate seeds for each call

    rng = np.random.default_rng(int(os.getpid()*time.time()))

    ins = [(int(n_samples/n_core)
            , cat_fn
            , void_fn
            , fp_fn
            , mem_lim/n_core
            , int(rng.uniform(0, int(time.time())
                              )
                 )
            ) for i in range(n_core)]
    
    with Pool(n_core) as p:
        simulated_data = p.starmap(monte_carlo, ins)
    # print(len(simulated_data))
    # print(len(simulated_data[0]))
    # print(len(simulated_data[0][0]))
    # master_list  = np.array([])
    master_list = pd.DataFrame()
    for dat in simulated_data:
        master_list = pd.concat((master_list, dat[0])).copy()
    print(f"Completed with exit codes: {[ dat[1] for dat in simulated_data]}")
    print("2 is what you want. 1 means you hit a memory cap")
    # Save the generated data
    save_dat(cat_fn, master_list)