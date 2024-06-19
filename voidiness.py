"""This holds the voidiness analysis seperate from the jupyter notebooks"""
from custom_functions import *
import pandas as pd

CEL_DATA_FN = "cel_obj_table.xlsx"
# Strictly Sticking to voidiness analysis we only need  positional data 'RAJ2000', 'DEJ2000', 'Redshift'

VOIDS_DATA_FN = "processed_voids.xlsx"
# Needs columns 'ID', 'RAdeg', 'DEdeg', 'z', 'Reff', 'cmvd_Mpc', 'Reff_Mpc', 'r_ang_deg'
# - ID: Void ID from catalog
# - RAdeg: Right Ascension in degrees in J2000 Frame
# - DEDeg: Declination "                              "
# - z: Redshift
# - cmvd_Mpc: Comoving distance to void in megaparsecs
# - Reff_Mpc: Effective on-sky effective radius " "
# - r_ang_deg: Void angular radius in degrees


# One of the first attempts. Keeping it for posterity
vhes = {}

bad_ints = {'v_idx': [],
            's_idx': []}
def pre_voidy_calc(voids, cel_obj):
    for v_idx in voids.index:
        temp_vhe = cel_obj.copy() # Save a fresh copy of work vhe

        # Grab void data
        void_ra, void_de, = voids.loc[v_idx,['RAdeg', 'DEdeg']]
        r_ang_deg, z, v_cmvd, v_r_mpc = voids.loc[v_idx,['r_ang_deg', 'z', 'cmvd_Mpc', 'Reff_Mpc']]

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
        else:
            continue


        for grs_idx in s_idx:
            ra, de, s_cmvd = cel_obj.loc[grs_idx, ['RAdeg', 'DEdeg', 'cmvd_Mpc']]
            singular_vhe_skycoord = SkyCoord(ra * u.deg, de*u.deg)

            s_v_dist = temp_void_coord.separation(singular_vhe_skycoord) # source-void distnace

            if s_v_dist.deg < r_ang_deg:
                # Last check to ensure sources are inside voids
            
                void_int, Cv_i, bad_int = calulate_voidy_int(void_ra, void_de, v_cmvd, v_r_mpc,
                                            ra, de, s_cmvd,
                                            s_v_dist)

                if bad_int:
                    bad_ints['v_idx'].append(v_idx)
                    bad_ints['s_idx'].append(grs_idx)
                data = vhes.setdefault(grs_idx, {
                                            'void_idx': [],
                                            'Cv': [],
                                            'intervals': []
                                        })    
                        

                data['void_idx'].append(v_idx)
                data['Cv'].append(Cv_i)
                data['intervals'].append(void_int)





def main():
    voids = pd.read_excel(VOIDS_DATA_FN)
    cel_obj = pd.read_excel(CEL_DATA_FN)
    pre_voidy_calc(voids,cel_obj)

if __name__ == "__main__":
    main()