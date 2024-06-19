"""This holds the voidiness analysis seperate from the jupyter notebooks"""
from custom_functions import *
import pandas as pd


# Strictly Sticking to voidiness analysis we only need  positional data 'RAJ2000', 'DEJ2000', 'Redshift'
CEL_DATA_FN = "cel_obj_table.xlsx"

# Needs columns 'ID', 'RAdeg', 'DEdeg', 'z', 'Reff', 'cmvd_Mpc', 'Reff_Mpc', 'r_ang_deg'
# - ID: Void ID from catalog
# - RAdeg: Right Ascension in degrees in J2000 Frame
# - DEDeg: Declination "                              "
# - z: Redshift
# - cmvd_Mpc: Comoving distance to void in megaparsecs
# - Reff_Mpc: Effective on-sky effective radius " "
# - r_ang_deg: Void angular radius in degrees

VOIDS_DATA_FN = "processed_voids.xlsx"



def main():
    voids = pd.read_excel(VOIDS_DATA_FN)
    cel_obj = pd.read_excel(CEL_DATA_FN)
    pass

if __name__ == "__main__":
    main()