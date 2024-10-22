import plotly.express as px
import pandas as pd
from custom_functions import make_astro_vec2
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
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

def plot_with_plotly(voids):

    # Removing lrg voids from sutter.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # lrgbright_idx = voids.index[voids['label'] == "2lrgbright" ]
    # lrgdim_idx = voids.index[voids['label'] == "2lrgdim"]
    # lrg_idx = lrgbright_idx.append(lrgdim_idx)

    # voids = voids[~voids['label'].isin(["2lrgbright", "2lrgdim"])]

    # voids = voids.iloc[lrg_idx]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mao_idx = voids.index[voids['label'] == "Mao"]
    mao_voids = voids.loc[mao_idx]

    sutter_idx = voids.index[~(voids['label'] == "Mao")]
    sut_voids = voids.loc[sutter_idx]

    # Add a column 'Dataset' to distinguish between sut_voids and mao_voids
    sut_voids['Dataset'] = 'Sut Voids'
    mao_voids['Dataset'] = 'Mao Voids'
    # Concatenate sut_voids and mao_voids into a single DataFrame
    combined_voids = pd.concat([sut_voids, mao_voids])

    #
    # Plots them together different colors
    #
    fig = px.scatter_3d(combined_voids, x='xMpc', y='yMpc', z='zMpc',
                    size='Reff_Mpc', color='Dataset',
                    opacity=0.7, title='Combined Voids')
    fig.show()

    #
    # Plots them together same color.
    #
    fig = px.scatter_3d(voids, x='xMpc', y='yMpc', z='zMpc',
                size='Reff_Mpc', 
                opacity=0.7)
    fig.show()


    
    # Plots them seperately
    
    fig1 = px.scatter_3d(sut_voids, x='xMpc', y='yMpc', z='zMpc',
                     size='Reff_Mpc', color_discrete_sequence=['blue'],
                     opacity=0.7, title='Sut Voids')

    # Plotting mao_voids in red
    fig2 = px.scatter_3d(mao_voids, x='xMpc', y='yMpc', z='zMpc',
                        size='Reff_Mpc', color_discrete_sequence=['red'],
                        opacity=0.7, title='Mao Voids')
    fig1.show()
    fig2.show()



def plot_matplotlib_dots():
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(voids.xMpc, voids.yMpc, voids.zMpc, alpha=0.2, s=voids.Reff_Mpc)
    plt.show()



def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:7j, 0:np.pi:3j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)


def plot_spheres_matplotlib():
    x, y, z = voids.xMpc, voids.yMpc, voids.zMpc
    r = voids.Reff_Mpc


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # draw a sphere for each data point
    for (xi,yi,zi,ri) in zip(x,y,z,r):
        (xs,ys,zs) = drawSphere(xi,yi,zi,ri)
        ax.plot_surface(xs, ys, zs, color=(0.6,0.8,1, 0.2))

    ax.set_aspect('equal')
    plt.show()

if __name__ == '__main__':
    # plot_spheres_matplotlib()
    plot_with_plotly(voids)
 