"""To determine a 2D footprint of the cosmic voids from the Sutter Catalog."""
from descartes import PolygonPatch # comes with bug. See NOTE below.
from matplotlib import pyplot as plt
import pandas as pd
import alphashape # used to determine the outer perimeter of a set of points


voids = pd.read_excel('exported_dataFrames/voids.xlsx')


# Generate tuples of coordinate pairs and save to a list becase its what alpha shape wants
v_coords = []
for ra, de in zip(voids.RAdeg, voids.DEdeg):
    v_coords.append((ra,de))

# Generate the alphashape
alpha = 0.20 # tightness of the shape. Higher is tighter. To high and data points are lost
a_shape = alphashape.alphashape(v_coords, alpha=alpha)
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
plt.scatter(ra.tolist(), de.tolist(), marker='.', color= 'red')
plt.show()