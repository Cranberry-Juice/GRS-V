"""This holds the voidiness analysis seperate from the jupyter notebooks"""
from custom_functions import *
import pandas as pd

CEL_DATA_FN = "cel_obj_table.xlsx"
# Strictly Sticking to voidiness analysis we only need  positional data 'RAJ2000', 'DEJ2000', 'Redshift'

VOIDS_DATA_FN = "processed_voids.xlsx"
# TODO: Double check the necessary columns.
# Needs columns 'ID', 'RAdeg', 'DEdeg', 'Reff', 'cmvd_Mpc', 'Reff_Mpc', 'r_ang_deg'
# - ID: Void ID from catalog
# - RAdeg: Right Ascension in degrees in J2000 Frame
# - DEDeg: Declination "                              "
# - z: Redshift
# - cmvd_Mpc: Comoving distance to void in megaparsecs
# - Reff_Mpc: Effective on-sky effective radius " "
# - r_ang_deg: Void angular radius in degrees


def pre_voidy_calc(voids, cel_obj):

    vhes = {}
    """Returns dict with cel_obj index as key entries. Each entry is associated
     with another dictionary with the intersecting void indices, the chord length of line of sight through the void and the mathematical interval of fraction of the LoS through void as entries. 

     NOTE: Indices refers to the row index of the data point in the pandas dataframe
    
    # 24 refers to the celestial object in row 24 of the inputted pandas dataframe
    EX. {24: {void_idx: [1,3,5,6], # Void indices 
              Cv: [0.1, 2, 0.2, 0.15], # Chord length of LoS and Void intersection Mpc
              interval: [ [0.62, 0.66], # LoS enters at 0.62 and exits at 0.66 
                                        along Los
                          [0.54, 0.55],
                          [0.44, 0.46],
                          [0.53, 0.57]
                        ]
              }
        }
                           # list of interval.Interval objects representing 
                           # entry and exit as fractions along LoS}"""
    # NOTE: Code originally written for gamma ray sources but is now generalized 
    # for any celestiabl object. Source in this function refers to celestial
    # object in question, ei galaxies, stars. etc.
    for v_idx in voids.index:
        temp_vhe = cel_obj.copy() # Save a fresh copy of celestial objects table

        # Grab void data
        void_ra, void_de, = voids.loc[v_idx,['RAdeg', 'DEdeg']]
        r_ang_deg, v_cmvd, v_r_mpc = voids.loc[v_idx,['r_ang_deg', 'cmvd_Mpc', 'Reff_Mpc']]

        temp_void_coord = SkyCoord(void_ra * u.deg, void_de *u.deg)
        vhe_coords =  SkyCoord(temp_vhe.RAdeg.values * u.deg, temp_vhe.DEdeg.values * u.deg)

        dist = temp_void_coord.separation(vhe_coords).deg

        radius_mask = dist  < r_ang_deg

        # If any grs are within radius of void
        if any(radius_mask):
            # Get indices of sources
            s_idx  = temp_vhe.index[radius_mask]
        else:
            continue # No grs within this void

        # Filter by having them be at least inside the void
        behind_mask = temp_vhe.loc[s_idx, 'cmvd_Mpc'] > (v_cmvd  - v_r_mpc)

        if any(behind_mask):
            s_idx = behind_mask.index[behind_mask]



        for grs_idx in s_idx:
            ra, de, s_cmvd = cel_obj.loc[grs_idx, ['RAdeg', 'DEdeg', 'cmvd_Mpc']]
            singular_vhe_skycoord = SkyCoord(ra * u.deg, de*u.deg)

            s_v_dist = temp_void_coord.separation(singular_vhe_skycoord) # source-void distnace

            if s_v_dist.deg < r_ang_deg:
                # Last check to ensure sources are inside voids
                void_int, Cv_i, bad_int = calulate_voidy_int(void_ra, void_de, v_cmvd, v_r_mpc,
                                            ra, de, s_cmvd,
                                            s_v_dist)
                
                data = vhes.setdefault(grs_idx, {
                                            'void_idx': [],
                                            'Cv': [],
                                            'intervals': []
                                        })
                
                data['void_idx'].append(v_idx)
                data['Cv'].append(Cv_i)
                data['intervals'].append(void_int)
    return vhes


def calc_master_voidiness(int_dict, cel_obj):
    # Add new column to save voidiness values for each source
    cel_obj = cel_obj.assign(Voidiness=np.zeros(len(cel_obj)))
    for idx in list(cel_obj.index):
        if idx in int_dict.keys():

            # Interval based voidiness calculation
            ints = int_dict[idx]['intervals']
            union = take_union(ints)
            cel_obj.at[idx, 'Voidiness'] = calc_voidiness(union)

        else:
            cel_obj.at[idx, 'Voidiness']  = 0

    return cel_obj


def voidy_analysis(voids_data_fn, cel_data_fn, indexed = True):
    """
    Returns pandas dataframe of cel_data table with Voidiness column appeneded. 
    Indexed means if the first column contains the indices from earlier pandas dataframe. True if first column is the indices from past dataframes"""
    voids = pd.read_excel(voids_data_fn)
    if indexed:
        cel_obj = pd.read_excel(cel_data_fn, index_col=0)
    else:
        cel_obj = pd.read_excel(cel_data_fn)
    
    intersect_data = pre_voidy_calc(voids, cel_obj)
    return calc_master_voidiness(intersect_data, cel_obj)


def main():
    results = voidy_analysis(VOIDS_DATA_FN, CEL_DATA_FN)
    plot_hist_watermark([results.Voidiness], bins=10,
          xlabel='Voidiness [$D_{\\rm Void}/D_{\\rm Total}$]',
          ylabel='Normalized Fraction',
          title='Voidiness Histogram',
          watermark_text='Preliminary',
          density=True)


if __name__ == "__main__":
    main()