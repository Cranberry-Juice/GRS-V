import plotly.express as px
import pandas as pd
from custom_functions import make_astro_vec2
import matplotlib.pyplot as plt


# Read in void data
voids = pd.read_excel('exported_dataFrames/voids.xlsx')

# Scale the void moving distance to be between 0 and 1
unit = max(voids.cmvd_Mpc)
voids.cmvd_Mpc = voids.cmvd_Mpc/unit
# Scale the radii as well
voids.Reff_Mpc = voids.Reff_Mpc/unit

# Convert to cartesian
xyz = make_astro_vec2(voids.cmvd_Mpc, voids.RAdeg, voids.DEdeg) # converts to cartesian
voids['xMpc'] = xyz[0]
voids['yMpc'] = xyz[1]
voids['zMpc'] = xyz[2]

# fig = px.scatter_3d(voids, x='xMpc', y='yMpc', z='zMpc',
#               size='Reff_Mpc', 
#               opacity=0.7)
# fig.show()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(voids.xMpc, voids.yMpc, voids.zMpc, alpha=0.2, s=voids.Reff_Mpc)
plt.show()