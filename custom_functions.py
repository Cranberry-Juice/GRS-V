"""The purpose of this script is to hold functions used to generate random 
SkyCoordsin the galactic plane and measure the distances between all possible sources"""

import numpy as np
from matplotlib import pyplot as plt
# import matplotlib.pyplot as plt
from itertools import chain
from astropy.coordinates import SkyCoord
from astropy import units as u
import math
from collections import defaultdict
import intervals as I # Used to do interval union math
import random
import pandas as pd
# import portion as I

def z_mirror(qsrs, binning, count_to_delete=None):
    idx_buckets = []
    # count_to_delete = [2.000e+00, 1.380e+02, 2.680e+02, 4.680e+02, 7.240e+02, 1.008e+03,
    #    1.307e+03, 1.811e+03, 2.211e+03, 2.898e+03, 3.654e+03, 3.902e+03]
    # count_to_delete = [   0.,  137.,  266.,  467.,  723., 1006., 1304., 1809., 2211.,
    #    2897., 3653., 3901.]
    # count_to_delete = [1.000e+00, 5.390e+02, 1.596e+03, 2.910e+03, 5.003e+03, 7.447e+03]
    if count_to_delete is None:
        count_to_delete = [   0., 1085., 1863., 2202., 3074., 3408.]
    # binning = np.linspace(0.4, 0.7, 7)
    for i in range(len(binning)-1):
        if i == len(binning)-2:
            too_high = qsrs.z > binning[i+1]    
        else:
            too_high = qsrs.z >= binning[i+1]
        too_low = qsrs.z < binning[i]
        bin_mask = ~too_high & ~too_low
        temp = qsrs[bin_mask]
        candidate_idx = temp.index.to_list()
        popped_idx = random.sample(candidate_idx, int(count_to_delete[i]))
        idx_buckets.append(popped_idx)
        qsrs.drop(popped_idx, axis='index', inplace=True)
    return (qsrs, idx_buckets)

# Get indices of simulated voidiness values that fall within specific redshift
def get_sim_idx(sim, real_df, zbin):
    mask = (real_df.z < zbin[1]) & (real_df.z>=zbin[0])
    idx_mask  = [i for i, val in enumerate(mask) if val]
    len_dat = len(real_df)
    l_sim = len(sim)
    n_samples = l_sim/len_dat

    assert n_samples%1 == 0, "Length of simulated datapoints is not integer multiple of length of data"

    n_samples = int(n_samples)
    # for i in range(n_samples):
    #     for idx in idx_mask:
    #         idx + len_dat*i
    return [idx + len_dat*i for i in range(n_samples) for idx in idx_mask]

    
    
def z_mask(df, zbin):
    return (df.z < zbin[1]) & (df.z >= zbin[0])

def calc_voidiness(union):
    voidiness = 0
    for voidichord in union:
        if ~voidichord.is_empty(): # This protects against empty intervals
            try:
                voidiness += voidichord.upper - voidichord.lower
            except:
                continue # get sidestepped 
            # Ok but actually, some intervals are atomic and have a neg inf and
            # pos inf upper bound. This is a deeper problem. 
            # Some ints, their lower bound is higher than the upper. This needs to
            # be investigated
            # TODO: investigate emtpy intervals.
            # I am back from the future. I checked it.
            # The empty intervals are actually properly detecting edge cases,
            # where the source lies just barely in front of the void, but passes
            # all tests as if it were inside or behind.
            # They will be kept in the code to detect and discard bad edge cases. 
    return voidiness

def take_union(intervals):
    """Returns the union of all the list inteverals."""
    union = I.open(0, 0)
    for thing in intervals:
        if thing.is_empty():
            # raise TypeError("interval is atomic") # it might not actually be a bug. This might be properly detecting edge cases
            # print(thing)
            continue
        else: # This protects against empty intervals
            union = union | thing
            # print(thing)
    return union

def w1(a, b):
    """Projecton of a onto b"""
    return (np.dot(a, b)/ mag(b)**2) * b

def cart_comp(r, th, ph):
    x = r * np.sin(th) * np.cos(ph)
    y = r * np.sin(th) * np.sin(ph)
    z = r * np.cos(th)

    return np.array([x, y, z])
def calulate_voidy_int(v_ra, v_de, v_cmvd, v_r_mpc, s_ra, s_de, s_cmvd, s_v_sep):

    # Collat
    v_th = 90 - v_de
    s_th = 90 - s_de

    v_th *= (math.pi)/180
    s_th *= (math.pi)/180
    

    
    v_ra *= math.pi/180
    s_ra *= math.pi/180


    # This was the bug
    rv = v_cmvd * np.array([1, v_th, np.sin(v_th)*v_ra])
    # rv += v_cmvd/mag(rv) # This should make it so i stop getting huge numbers

    rv_cart = cart_comp(v_cmvd, v_th, v_ra)
    rs_cart = cart_comp(s_cmvd, s_th, s_ra)

    # Calculating intersection and intersectoin intervals
    w1_vec = w1(rv_cart, rs_cart)
    # d_vec = rv_cart - w1_vec
    
    # Cv_i = 2 *np.sqrt(v_r_mpc**2 - mag(d_vec)**2)
    little_d = s_v_sep.rad * v_cmvd # Arclength approximation

    ##
    Cv_i = other_voidy_calc(little_d, v_r_mpc) # Cv_i no longer nanny
    ##

    ## Normalization working properly
    w2 = Cv_i/2 * (rs_cart / mag(rs_cart))
    ## 

    T1 = w1_vec - w2

    T2 = T1 + 2 * w2


    n1 = mag(T1)/s_cmvd
    n2 = mag(T2)/s_cmvd

    if n2 > 1:
        n2 = 1

    interval = I.closed(n1, n2)
    if interval.is_empty():
        bad_int = True
    else:
        bad_int = False
    return I.closed(n1, n2), Cv_i, bad_int

def get_void_vec(v_idx, void_table):
    v_RA = void_table.loc[v_idx]['RAdeg']
    v_DE = void_table.loc[v_idx]['DEdeg']
    v_cmvd = void_table.loc[v_idx]['cmvd_Mpc']
    return make_astro_vec2(v_cmvd, v_RA, v_DE)

def get_mask(vhe_table, RA_lo=0, RA_hi=360, DE_lo=-90, DE_hi=90):
    """Custom function to get mask of data within ranges specified"""
    RA_msk1 = vhe_table.RAdeg > RA_lo
    RA_msk2 = vhe_table.RAdeg < RA_hi
    DE_msk1 = vhe_table.DEdeg > DE_lo
    DE_msk2 = vhe_table.DEdeg < DE_hi

    RA_mask = np.logical_and(RA_msk1, RA_msk2)
    DE_mask = np.logical_and(DE_msk1, DE_msk2)
    return np.logical_and(RA_mask, DE_mask)

def get_ands(masks):
    """Get bitwise and of all masks passed in."""
    prev_mask = masks[0]
    for i in range(1, len(masks)):
        master_mask = np.logical_and(prev_mask, masks[i])
        prev_mask = master_mask
    return prev_mask

def merge_dicts(d1, d2):
    dd = defaultdict(list)
    for d in (d1, d2): # you can list as many input dicts as you want here
        for key, value in d.items():
            dd[key].append(value)
    return dict(dd)

def c(ang):
    return np.cos(ang) # Makes code reasonable
def s(ang):
    return np.sin(ang)


def transform(v, th, ph, to):
    """
    v: vector in cartesian
    th: (rad) colatitidue. measured from z
    ph: (rad) longitude. measured from x
    to: target coordinate system ("cart" or "sphere") """
    rot_mat = np.array([[s(th) * c(ph),  c(th) * c(ph), -s(ph)],
                        [s(th) * s(ph),  c(th) * s(ph),  c(ph)],
                        [c(th)        , -s(th)        , 0    ]])
    
    if to == "cart":
        return np.matmul(rot_mat, v)
    elif to == "sphere":
        return np.matmul(rot_mat.transpose, v)
    else:
        raise ValueError("Invalid 'to' value. It must be either 'cart' or 'sphere'.")


def vector_transform(v, th, ph, fromto):
    v1 = v[0]
    v2 = v[1]
    v3 = v[2]
    if fromto == 'cart':
        va = v1 * s(th) * c(ph) + v2 * c(th) * c(ph) + v3 * -s(ph) 
        vb = v1 * s(th) * s(ph) + v2 * c(th) * s(ph) + v3 * c(ph)
        vc = v1 * c(th)         + v2 * -s(th)        + v3 * 0
    if fromto == 'sphere':
        va = v1 * s(th) * c(ph) + v2 * s(th) * s(ph) + v3 *  c(th) 
        vb = v1 * c(th) * c(ph) + v2 * c(th) * s(ph) + v3 * -s(ph)
        vc = v1 * -s(ph)        + v2 * c(th)         + v3 *  0

    return np.array([va, vb, vc])

def other_voidy_calc(little_d_mp , void_reff_Mpc):
    return 2* np.sqrt(void_reff_Mpc**2 - little_d_mp**2)

def calc_void_chord(rv, rs, void_r, ang_dist):
    """This is my contribution :^) """
    # dot1 = np.dot(rs,rv)
    # dot2 = np.dot(rs,rs)
    # w1 = (dot1/dot2) * rs 
    # d = np.sqrt(mag(rv)**2 - mag(w1)**2)
    # d = calc_d(rs, rv)
    d = np.sin(ang_dist*(math.pi/180)) * mag(rv)
    chord = 2 * np.sqrt(void_r**2 - d**2)
    # if d > void_r:
    #     print('imaginaries')
    return chord

def calc_void_chord2(rv, rs, void_r):
    """This is my contribution :^) """
    dot1 = np.dot(rs,rv)
    dot2 = np.dot(rs,rs)
    w1 = (dot1/dot2) * rs 
    d = np.sqrt(mag(rv)**2 - mag(w1)**2)
    # d = calc_d(rs, rv)
    # d = np.sin(ang_dist*(math.pi/180)) * mag(rv)
    chord = 2 * np.sqrt(void_r**2 - d**2)
    # if d > void_r:
    #     print('imaginaries')
    return chord
    
def voidy_inside(rs, rv, void_reff_mpc):
    mag_rs = mag(rs)
    mag_rv = mag(rv)
    th = angle_from_dot(rs,rv)
    phi = np.arcsin((np.sin(th) * mag_rs)/void_reff_mpc)

    rw = rv - rs
    mag_rw = mag(rw)
    beta = np.arccos((np.dot(rv,rw)/(mag_rv*mag_rw)))

    psi = beta - phi
    alf = np.arccos(np.dot(rw,rs)/mag_rw*mag_rs)
    cv = (void_reff_mpc * np.sin(psi))/np.sin(alf)
    return cv

def angle_from_dot(v1, v2):
    return np.arccos(np.dot(v1, v2)/(mag(v1)*mag(v2)))

def calc_d1(rs, rv):
    th = angle_from_dot(rs, rv) # rad
    return np.tan(th) * mag(rv)

def calc_d2(angle_rad, rv):
    return np.sin(angle_rad) * mag(rv)

def calc_d_inside(rs, rv):
    rw = rs - rv
    return mag(rw)

def make_astro_vec(cmvd, lon, lat):
    """
    Takes cmvd, lon and lat in degs returns spherical vector"""
    th = to_colat_rad(lat)
    ph = lon* (math.pi/180)
    return make_sph_vec(cmvd, th, ph)

def make_astro_vec2(cmvd, lon, lat):
    """
    Takes cmvd, lon and lat in degs returns cart nparray vector"""
    r = cmvd
    th = to_colat_rad(lat)
    ph = lon * (math.pi/180)
    x = r * np.sin(th) * np.cos(ph)
    y = r * np.sin(th) * np.sin(ph)
    z = r * np.cos(th)
    return np.array([x, y, z])

def make_sph_vec(r, th, ph):
    return np.array([r, r*th, r*np.sin(th)*ph])

def to_colat_rad(lat_deg):
    """Returns latitude in radians"""
    return (90 - lat_deg) * (math.pi/180)  
def add_to_dict(dict, key, value):
    try:
        dict[key] += [value]
    except:
        dict[key] = [value]
    return dict

def verify_combo_count(array_1d,n,k=2):
    """Raises error if n choose k number of combinations does not match
    array length and compares against length of array_1d"""
    n_combinations = math.factorial(n)/(math.factorial(k) * math.factorial(n - k))
    if n_combinations != len(array_1d):
        raise Exception(f"Error in producing combinations.\n Array length:{n}\n Num Combos: {n_combinations}")
    
    
def plot_gal_sources(ang_seps_deg):
    # ang_seps_deg = ang_seps_gal * ((180 * u.deg)/math.pi) This is a bad line
    plt.hist(ang_seps_deg, bins=100, histtype="step") # 100 bins chosen arbitrarily
    plt.xlabel(f'Distance {u.deg}')
    plt.ylabel('Count')
    plt.title('Angular Distance Between Sources (Galactic frame)')
    plt.grid()
    plt.tight_layout()
    plt.show()


def plot_hist(data, bins=180, 
              histtype="step", 
              grid=True,
              density=True, 
              stacked=False, 
              label=None,
              **kwargs):
    """
    Create a histogram plot using matplotlib.

    Keyword Arguments:
        xlabel (str, optional): Label for the x-axis.
        ylabel (str, optional): Label for the y-axis.
        title (str, optional): Title for the plot.
    """
    plt.hist(data, 
             bins=bins, 
             histtype=histtype, 
             density=density, 
             stacked=stacked,
             label = label)
    if 'xlabel' in kwargs:
        plt.xlabel(kwargs['xlabel'])
    if 'ylabel' in kwargs:
        plt.ylabel(kwargs['ylabel'])
    if 'title' in kwargs:
        plt.title(kwargs['title'])
    if grid:
        plt.grid()
    if label != None:
        plt.legend()
    plt.tight_layout()
    
    plt.show()


# import matplotlib.pyplot as plt

def plot_hist_watermark(data, bins=180, 
              histtype="step", 
              grid=True,
              density=True, 
              stacked=False, 
              label=None,
              watermark_text=None,  # Add a watermark text parameter
              **kwargs):
    plt.hist(data, 
             bins=bins, 
             histtype=histtype, 
             density=density, 
             stacked=stacked,
             label=label)
    
    if 'xlabel' in kwargs:
        plt.xlabel(kwargs['xlabel'])
    if 'ylabel' in kwargs:
        plt.ylabel(kwargs['ylabel'])
    if 'title' in kwargs:
        plt.title(kwargs['title'])
    if grid:
        plt.grid()
    if label is not None:
        plt.legend()
    
    if watermark_text is not None:
        plt.text(0.5, 0.5, watermark_text,  # Adjust the position as needed
                 fontsize=30, color='gray', alpha=0.5,
                 ha='center', va='center', transform=plt.gca().transAxes)
    
    plt.tight_layout()
    plt.show()





def sample_and_measure(n_coords):
    """Samples Random Coordinates and measures their angular distances,
    RE
    n_coords: Number of coordidate pairs to generate"""
    rand_coords_gal = generate_rand_Coords_static(n_coords)
    ang_seps_gal = dist_w_skycoord(rand_coords_gal)
    return ang_seps_gal

def len_gen_coords(n):
    return len(generate_rand_Coords_dynamic(n))

# global COUNT

# COUNT = 0
def generate_rand_Coords_static(n, mask=10, seeded=False):

    """Generates random galactic coordinates accounting for galactic plane bias with a plus_minus "mask"
     mask: IN DEGREES
     n: number of coordinate pairs to produce
     Returns: SkyCoord Array with n coordinates in the galactic frame """
    
    # pylint: disable=trialing-whitespace
    l_coords, b_coords = generate_rand_l_and_b(n + 211,seeded)
    rand_sources = SkyCoord(l=l_coords, b=b_coords, frame = "galactic")

    # Apply mask
    if mask != 0:
        # Find fermi data with |b| < 10
        mask_array = abs(rand_sources.b.value) <= mask

        # Keep the sources with |b| > 10
        rand_sources = rand_sources[~mask_array] # Drops about 50 sources
        # b_coords = b_coords[~mask_array]
        # l_coords = l_coords[~mask_array]

        # # Keep generating and masking sources until we match the intended number
        # # of data points
        # n_rand_sources = len(rand_sources)


        # add_b_coords = [] * u.deg # Avoids bug when appending deg with this list
        # add_l_coords = [] * u.deg
        # Generate more coords to make up for lost coords
        # while n_rand_sources < n: 
        #     n_missing = (n - n_rand_sources) * 3 # This should help speed up loop

        #     temp_l_coords, temp_b_coords = generate_rand_l_and_b(n_missing)

        #         # Find fermi data with |b| < 10
        #     temp_mask_array = abs(temp_b_coords.value) <= mask

        #     # Save additional sources with |b| > 10
        #     add_b_coords = np.append(add_b_coords,temp_b_coords[~temp_mask_array]) 
        #     add_l_coords = np.append(add_l_coords,temp_l_coords[~temp_mask_array])

        #     # Update new length of array
        #     n_rand_sources = len(rand_sources) + len(add_b_coords)

        # if n_rand_sources > n:
        #     n_surp = n_rand_sources - n
        #     add_b_coords = add_b_coords[:-n_surp]
        #     add_l_coords = add_l_coords[:-n_surp]
        #     # print("There were surplus coordinates")
        # # else:
        #     # print("No Surplus Coordiantes") 
        
        # # We should have the exact amount needed by now. 
        # b_coords = np.append(b_coords, add_b_coords)
        # l_coords = np.append(l_coords, add_l_coords)
        # rand_sources = SkyCoord(l=l_coords, b=b_coords,frame = "galactic")
    

    # i am possibly being too paranoid with this one
    # if len(rand_sources) != n:
    #     raise Exception("Code fail to generate n number of coordinates with\n\
    #                     required mask.")
    # print(COUNT)
    # COUNT +=1
    return  rand_sources

def generate_rand_Coords_dynamic(n, mask=10, seeded=False):

    """Generates random galactic coordinates accounting for galactic plane bias with a plus_minus "mask"
     mask: IN DEGREES
     n: number of coordinate pairs to produce
     Returns: SkyCoord Array with n coordinates in the galactic frame """
    
    # pylint: disable=trialing-whitespace
    l_coords, b_coords = generate_rand_l_and_b(n,seeded)
    rand_sources = SkyCoord(l=l_coords, b=b_coords, frame = "galactic")

    # Apply mask
    if mask != 0:
        # Find fermi data with |b| < 10
        mask_array = abs(rand_sources.b.value) <= mask

        # Keep the sources with |b| > 10
        rand_sources = rand_sources[~mask_array] # Drops about 50 sources
        b_coords = b_coords[~mask_array]
        l_coords = l_coords[~mask_array]

        # Keep generating and masking sources until we match the intended number
        # of data points
        n_rand_sources = len(rand_sources)


        add_b_coords = [] * u.deg # Avoids bug when appending deg with this list
        add_l_coords = [] * u.deg
        # Generate more coords to make up for lost coords
        while n_rand_sources < n: 
            n_missing = (n - n_rand_sources) * 3 # This should help speed up loop

            temp_l_coords, temp_b_coords = generate_rand_l_and_b(n_missing)

                # Find fermi data with |b| < 10
            temp_mask_array = abs(temp_b_coords.value) <= mask

            # Save additional sources with |b| > 10
            add_b_coords = np.append(add_b_coords,temp_b_coords[~temp_mask_array]) 
            add_l_coords = np.append(add_l_coords,temp_l_coords[~temp_mask_array])

            # Update new length of array
            n_rand_sources = len(rand_sources) + len(add_b_coords)

        if n_rand_sources > n:
            n_surp = n_rand_sources - n
            add_b_coords = add_b_coords[:-n_surp]
            add_l_coords = add_l_coords[:-n_surp]
            # print("There were surplus coordinates")
        # else:
            # print("No Surplus Coordiantes") 
        
        # We should have the exact amount needed by now. 
        b_coords = np.append(b_coords, add_b_coords)
        l_coords = np.append(l_coords, add_l_coords)
        rand_sources = SkyCoord(l=l_coords, b=b_coords,frame = "galactic")
    


def generate_rand_l_and_b(n, seeded=False, seed = 567307250717):
    """
    Generates l and b galactic coordinates in degrees
    """
    if seeded:
        seed = seed # Using seed while debugging
    else:
        seed = None
    rng = np.random.default_rng(seed) # no seed
    theta = np.arccos(1 - 2 * rng.uniform(0, 1, n)) * (180/math.pi) * u.deg # COLAT
    b = 90* u.deg - theta
    l = rng.uniform(0, 360, n) * u.deg
    return (l, b)

def generate_rand_l_and_b_old(n, seeded=False, seed = 567307250717):
    """
    Generates l and b galactic coordinates in degrees
    """
    if seeded:
        seed = seed # Using seed while debugging
    else:
        seed = None
    rng = np.random.default_rng(seed) # no seed

    x = rng.uniform(-1,1+ pow(10,-8),n)
    y = rng.uniform(-1,1+ pow(10,-8),n)
    z = rng.uniform(-1,1+ pow(10,-8),n)

    # Quickly checking the xyz distributions are different and independent
    # plot_hist([x, y, z], 
    #           title = 'x y z Distributions', 
    #           label=['x', 'y', 'z'],
    #           bins = 50)

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r) * (180/math.pi)* u.deg # theta, spherical Co-latitude
    b = 90 * u.deg - theta # Galactic latitude
    l = np.sign(y) * np.arccos(x/np.sqrt(x**2 + y**2)) * (180/math.pi) * u.deg # phi
    return (l, b)

def generate_rand_Coords_old(n, mask=10, seeded=False):
    """Generates random dgalactic coordinates accounting for galactic plane bias with a plus_minus "mask"
     n: number of coordinate pairs to produce
     Returns: SkyCoord Array with n coordinates in the galactic frame """
    
    # pylint: disable=trialing-whitespace
    if seeded:
        seed = 567307250717 # Using seed while debugging
    else:
        seed = None
    rng = np.random.default_rng(seed) # no seed

    n_total = n 

    plus_minus = 10 # ± 10 Degrees from the galactic plane. Kept for legacy purposes

    # plus_minus = plus_minus * (math.pi/180) # Already converted to rads later
    
    if n_total % 2 != 0: # n_total is odd
        # Randomly assign the extra source to the north or south hemisphere
        if rng.choice([True, False]):
            n_south = int(n_total//2  + 1) # South gets extra source
            n_north = int(n_total//2)
        else:
            n_south = int(n_total//2) 
            n_north = int(n_total//2 + 1) # North gets extra source
    # n_total is even
    else:               
        n_south = int(n_total/2)
        n_north = n_south

    # Check everything adds up
    if n_total != n_south + n_north:
        # pylint:  disable=broad-exception-raised
        raise ValueError("South and north numbers are not totalling to n_total")

    # # Generate random galactic coordinates in degrees


    # # Bounds on galactic latitude
    # north_low_bound = plus_minus # 10 Degs
    # north_upp_bound = 90 + pow(10, -8) # 90 Deg
    # # rng.uniform excludes upper bound but 90.000 000 00 is a valid coord
    # # This puts upper bound at 90.000 000 01 and should allow for 90 to be randomly picked

    # south_low_bound = -90 # -90 deg
    # south_upper_bound = -plus_minus # Degs
    # south_upp_bound = south_upper_bound + pow(10,-8) # -10 Deg


    # north_b = rng.uniform(north_low_bound, north_upp_bound, size=n_north) # 10 -> 90
    # south_b = rng.uniform(south_low_bound, south_upp_bound, size=n_south) # -90 -> -10

    # north_l = rng.uniform(0, 360, size=n_north)
    # south_l = rng.uniform(0, 360, size=n_south)
    
    # b_coords = np.append(north_b, south_b) * u.deg # DEG
    # l_coords = np.append(north_l, south_l) * u.deg # DEG

    b_coords = rng.uniform(-90, 90 + pow(10,-8), n_total) * u.deg
    l_coords = rng.uniform(0, 360, n_total) * u.deg
    rand_sources = SkyCoord(l=l_coords, b=b_coords,frame = "galactic")

    if mask != 0:
        # Find fermi data with |b| < 10
        mask_array = abs(rand_sources.b.value) <= mask

        # Keep the sources with |b| > 10
        rand_sources = rand_sources[~mask_array] # Drops about 50 sources
        b_coords = b_coords[~mask_array]
        l_coords = l_coords[~mask_array]

        # Keep generating and masking sources until we match the intended number
        # of data points
        n_rand_sources = len(rand_sources)
        add_b_coords = [] * u.deg # Avoids bug when appending deg with this list
        add_l_coords = [] * u.deg

        # Generate more coords to make up for lost coords
        while n_rand_sources < n_total: 
            n_missing = n_total - n_rand_sources

            temp_b_coords = rng.uniform(-90, 90 + pow(10,-8), n_missing) * u.deg
            temp_l_coords = rng.uniform(0, 360, n_missing) * u.deg
                # Find fermi data with |b| < 10
            temp_mask_array = abs(temp_b_coords.value) <= mask

            # Save additional sources with |b| > 10
            add_b_coords = np.append(add_b_coords,temp_b_coords[~temp_mask_array]) 
            add_l_coords = np.append(add_l_coords,temp_l_coords[~temp_mask_array])

            # Update new length of array
            n_rand_sources = len(rand_sources) + len(add_b_coords)

        # Remove surplus entries.
        # likely this is redundant code. 
        # i added an exception to cover this JUST in case.
        # if n_rand_sources > n_total:
        #     n_surp = n_rand_sources - n_total
        #     add_b_coords = add_b_coords[:-n_surp]
        #     add_l_coords = add_l_coords[:-n_surp]
        #     print("There were surplus coordinates")
        # else:
        #     print("Redundant")

        # We should have the exact amount needed by now. 
        b_coords = np.append(b_coords, add_b_coords)
        l_coords = np.append(l_coords, add_l_coords)
        rand_sources = SkyCoord(l=l_coords, b=b_coords,frame = "galactic")

    # i am possibly being too paranoid with this one
    if len(rand_sources) != n_total:
        raise Exception("Code fail to generate n number of coordinates with\n\
                        required mask.")
    return  rand_sources

def mag(vec):
    """ Returns magnitude of Vector"""
    return np.sqrt(sum(pow(element, 2) for element in vec))

def dist_vincenty(sky_coords):
    """
    Calculates the on-sky angular distance between all SkyCoord cordinates 
    using Vinencenty's formula. 

    Optionally takes in galactic frame input. Only accepts "gal" string.
        Otherwise it will use RA,DEC

    Returns: 1D array of angular seperations in radians 
    """
        
    vhe_coords = sky_coords
    # Originally code was meant for vhe_coords. Repurposed code into a function
    ang_seps_rad = []
    for i in range(1, len(vhe_coords)):
        frontCoords = vhe_coords[:-i]
        backCoords  = vhe_coords[i:]

        for j in range(len(frontCoords)):
            a = frontCoords[j] # Coordinate pair a
            b = backCoords[j] # Coordinate pair b
            # In spherical:
            # RA --> phi, 
            # DEC --> Theta  https://wikimedia.org/api/rest_v1/media/math/render/svg/b19bc4d661d9cf7cd09027568387e46ee452953f
            # From galactic:
            # l --> phi, 
            # b --> theta 
            r = 1 # All coords are put on surface of unit sphere
            phi_a = a.spherical.lon.rad
            phi_b = a.spherical.lat.rad

            theta_a = a.spherica.lon.rad
            theta_b = b.spherical.lat.rad

            # NOTE  on l latitude and DEC. Spherical coord theta is a 
            # co-latitude. Conversion is necessary. 

            # COMMENTED OUT FOR VINCENTY. Comment back in for vectors
            # PI = math.pi
            # theta_a = PI/2 - theta_a
            # theta_b = PI/2 - theta_b
            # A_vec = [r, r * theta_a, r * np.sin(theta_a) * phi_a ]
            # B_vec = [r, r * theta_b, r * np.sin(theta_b) * phi_b ]

            # This  might be the instability
            # accos method
            # temp_ang = np.arccos(np.dot(A_vec, B_vec) 
            #                      /(mag(A_vec) * mag(B_vec))) # rad 
            
            # atan method
            # temp_ang = np.arctan(mag(np.cross(A_vec, B_vec))/np.dot(A_vec, B_vec))

            # # Euclidian
            # d_phi = phi_a - phi_b
            # d_theta = theta_a = theta_b
            # temp_ang = np.sqrt(d_phi**2 + d_theta**2)

            # Vincenty Formula
            l1 = phi_a
            l2 = phi_b
            phi1 = theta_a
            phi2 = theta_b
            d_lat = l2 - l1
            t1 = np.cos(phi2) * np.sin(d_lat)
            t2 = np.cos(phi1) * np.sin(phi2) - \
                 np.sin(phi1) * np.cos(phi2) * np.cos(d_lat)
            opp = np.sqrt(t1**2 + t2**2)
            adj = np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(d_lat)
            temp_ang = np.arctan2(opp, adj)
            ang_seps_rad.append(temp_ang)
    return ang_seps_rad

def dist_w_skycoord(sky_coords):
    """Calculates the on-sky angular distance between all SkyCoord cordinates
     using SkyCoord.separation() method.
      Returns: 1D array of angular seperation values (unitless) """
        
    vhe_coords = sky_coords # Originally code was meant for vhe_coords. Repurposed code into a function
    ang_seps = []
    # USING ARRAYS. This method is significantly faster
    for i in range(1, len(vhe_coords)): # 0 index skipped. Gives dist. to self
        frontCoords = vhe_coords[:-i]
        backCoords  = vhe_coords[i:]
        temp_seps = frontCoords.separation(backCoords)
        temp_seps = temp_seps.value
        ang_seps.append(temp_seps)

    ang_seps = list(chain.from_iterable(ang_seps)) # Flattens nested lists
    '''
    The verification calculates factorials and is prolly just slowing down
    the code for no reason. It has been tested several times before.
    Code appears to properly create all combinations without double counting.
    '''
    # verify_combo_count(ang_seps, len(vhe_coords)) # if quiet we are good. 

# Going through 1 by one # Again, it makes no difference to histogram plot
    # for i in range(1, len(vhe_coords)):
    #     frontCoords = vhe_coords[:-i]
    #     backCoords  = vhe_coords[i:]
    #     for j in range(len(frontCoords)):
    #         a = frontCoords[j]
    #         b = backCoords[j]
    #         temp_sep = a.separation(b)
    #         ang_seps.append(temp_sep)

    return ang_seps

    