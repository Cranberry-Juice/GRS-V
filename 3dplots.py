import plotly.express as px
import pandas as pd


# Read in void data
voids = pd.read_excel('exported_dataFrames/voids.xlsx')

# Scale the void moving distance to be between 0 and 1
unit = max(voids.cmvd_Mpc)
voids.cmvd_Mpc = voids.cmvd_Mpc/unit
# Scale the radii as well
voids.Reff_Mpc = voids.Reff_Mpc/unit