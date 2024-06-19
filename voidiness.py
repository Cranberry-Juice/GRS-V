"""This holds the voidiness analysis seperate from the jupyter notebooks"""
from custom_functions import *
import pandas as pd

CEL_DATA_FN = "Cel_obj_data.xlsx"
VOIDS_DATA_FN = "processed_voids.xlsx"

def import_voids(fn):
    """
    Read void data from a text file into a pandas DataFrame.

    Parameters:
    fn (str): File path to the void data file.

    Returns:
    pandas.DataFrame: DataFrame containing void data with columns:
        - ID: Void ID from catalog
        - RAdeg: Right Ascension in degrees in J2000 Frame
        - DEDeg: Declination "                              "
        - z: Redshift
        - cmvd_Mpc: Comoving distance to void in megaparsecs
        - Reff_Mpc: Effective on-sky effective radius " "
        - r_ang_deg: Void angular radius in degrees
    """
    # Read void table data from text file
    voids = pd.read_table(fn, sep=",")  # Void data already parsed
    voids.drop(columns=voids.columns[0], axis=1, inplace=True) # Removes redundant row counter column
    return voids
    # NOTE: Unused columns: 'ID', 'RAdeg', 'DEdeg', 'z', 'Reff', 'cmvd_Mpc','Reff_Mpc', 'r_ang_deg'

def main():
    voids = import_voids(VOIDS_DATA_FN)
    pass

if __name__ == "__main__":
    main()