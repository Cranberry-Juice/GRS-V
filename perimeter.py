"""To determine a 2D footprint of the cosmic voids from the Sutter Catalog."""
from descartes import PolygonPatch # comes with bug. See NOTE below.
from matplotlib import pyplot as plt
import pandas as pd
import alphashape # used to determine the outer perimeter of a set of points
import numpy as np


voids = pd.read_excel('exported_dataFrames/voids.xlsx')


# Generate tuples of coordinate pairs and save to a list becase its what alpha shape wants
v_coords = []
for ra, de in zip(voids.RAdeg, voids.DEdeg):
    v_coords.append((ra,de))

# Generate the alphashape
alpha = 0.20 # tightness of the shape. Higher is tighter. To high and data points are lost

def variable_alpha(indices, r):
    top_alpha = 0.2 # alpha value for the top of the data points
    bottom_alpha = 0.4
    # filter by declination. Points with DE higher than 20 dec get top alpha
    # lower than 20 dec gets bottom alpha
    cutoff_de = 20 #deg

    return top_alpha if any(np.array(v_coords)[indices][:, 1] > cutoff_de) else bottom_alpha

a_shape = alphashape.alphashape(v_coords, variable_alpha)
fig, ax = plt.subplots()
ax.scatter(*zip(*v_coords), marker='.')


# Plot alpha shape
#
# NOTE !!!!!!! THIS REQUIRED ME TO MANUALLY EDIT THE DESCARTES PACKAGE IN MY ENVIRONMENT
# Link to fix
# https://stackoverflow.com/questions/75287534/indexerror-descartes-polygonpatch-wtih-shapely
# Otherwise this will not plot

ax.add_patch(PolygonPatch(a_shape, alpha=0.2))
ax.set_ylabel('Declination (Deg)')
ax.set_xlabel('Right Ascencion (Deg)')
ax.set_title('Cosmic Void Centers and Perimeter')

# Mark the void centers that define the boundary
ra, de = a_shape.exterior.coords.xy
ra = ra.tolist()
de = de.tolist()

# We have raw RA, DE but we need to know what voids they are to get radius information
maskra = voids.isin(ra).RAdeg
maskde = voids.isin(de).DEdeg

# Combine the mask. len(mask) < len(ra) by one. This is expected since the first
# ra, dec pair is also the last to define the polygon. 
# AKA fence post error, though not an error in this case.
mask = [x and y for x, y in zip(maskra, maskde)]

idx = voids[mask].index # indices of voids

# plot in a different color to confirm
ax.scatter(voids.loc[idx, "RAdeg"], voids.loc[idx, "DEdeg"], marker=".", color = 'red')
plt.show()
