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


# generate points at the circumference of the perimeter voids

n_samples = 10 # number of samples around the perimeter of the void
sampled_ra = np.zeros(len(idx * n_samples))
sampled_de = np.zeros(len(idx * n_samples))

thetas = np.linspace(0, 2*np.pi, n_samples)

oofra = np.array([])
oofde = np.array([])
for ix in idx:
    r_ang = voids.loc[ix,'r_ang_deg']
    x, y = voids.loc[ix, "RAdeg"], voids.loc[ix, "DEdeg"]

    ras = r_ang * np.cos(thetas) + x # + x for coordinate shift
    des = r_ang * np.sin(thetas) + y
    oofra = np.append(oofra, ras)
    oofde = np.append(oofde, des)


def tuple_and_list(x, y):
    """Take array of x and y's and return them as a list of (x,y) tuples"""
    llllll = []
    for xx, yy in zip(x, y):
        llllll.append((xx,yy))
    return llllll
# ax.scatter(oofra, oofde, marker='.')
# On second thought this might not be the best way to go about it but I will go
# through with it anyways to confirm or deny suspicions

coords = []
for raa, dee in zip(oofra, oofde):
    coords.append((raa, dee))
bigger_border = alphashape.alphashape(coords, variable_alpha)

ax.add_patch(PolygonPatch(bigger_border, alpha=0.2))

# Sample ALL THE VOIDS
n_samples = 60 # number of samples around the perimeter of the void
thetas = np.linspace(0, 2*np.pi, n_samples)

all_sample_ra = np.array([])
all_sample_de = np.array([])
for index in voids.index:
    r_ang = voids.loc[index,'r_ang_deg']
    x, y = voids.loc[index, "RAdeg"], voids.loc[index, "DEdeg"]

    ras = r_ang * np.cos(thetas) + x # + x for coordinate shift
    des = r_ang * np.sin(thetas) + y
    all_sample_ra = np.append(all_sample_ra, ras)
    all_sample_de = np.append(all_sample_de, des)

# ax.scatter(all_sample_ra, all_sample_de, marker='.')

all_coords = tuple_and_list(all_sample_ra, all_sample_de)

def variable_alpha2(indices, r):
    top_alpha = 0.4 # alpha value for the top of the data points
    bottom_alpha = 0.1
    # filter by declination. Points with DE higher than 20 dec get top alpha
    # lower than 20 dec gets bottom alpha
    upper_de = -4 #deg

    return top_alpha if any(np.array(all_coords)[indices][:,1]>upper_de) else bottom_alpha


all_perimeter = alphashape.alphashape(all_coords, alpha=0.11)

# Get set of points that defined the boundary
ra, de = all_perimeter.exterior.coords.xy
ra = ra.tolist()
de = de.tolist()

ax.scatter(ra,de, marker='.')


ax.add_patch(PolygonPatch(all_perimeter, alpha=0.2))
# plt.show()

#
# export the data points that define the boundary
#
defining_pts = {'RAdeg': ra,
                'DEdeg': de}

defining_pts = pd.DataFrame(defining_pts)
defining_pts.to_excel('exported_dataFrames/footprint_points.xlsx', index=False)
