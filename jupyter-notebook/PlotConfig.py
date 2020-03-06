#!/usr/bin/env python
# coding: utf-8

# In[329]:


import numpy as np
import h5py
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3


# In[32]:


f = h5py.File('./spaceConfiguration.h5','r')


# In[33]:


dset = f['links'][0]
k = f['coordinates'][0][()]
f.close()


# In[150]:
#
# Define arbitary dodecahedron and icosahedron
#

def initializeDodecahedron():
    phi = (1.+math.sqrt(5.))/2.
    alpha = math.asin(phi/math.sqrt(3.))-math.acos(phi/math.sqrt(phi+2.))
    beta = math.atan(2.*phi**2)
    theta = np.empty(20,dtype = np.float64)
    psi = []
    
    theta[0:5] = alpha
    theta[5:10] = beta
    theta[10:15] = math.pi - alpha
    theta[15:20] = math.pi - beta
    
    inc = np.arange(1./5.,10./5.,2./5.)*math.pi
    psi.append(inc)
    psi.append(inc)
    inc = np.arange(6./5.,3.,2./5.)*math.pi
    psi.append(inc)
    psi.append(inc)
    psi = np.asarray(psi,dtype = np.float64).flatten()
    return theta,psi

def initializeIcosahedron():
    gamma = math.pi/2. - math.atan(1./2.)
    
    theta = np.zeros(12, dtype = np.float64)
    theta[1:6] = gamma
    theta[6] = math.pi
    theta[7:12] = math.pi - gamma
    
    psi = []
    inc = np.arange(0,9./5.,2./5.)*math.pi
    psi.append(inc)
    inc = np.arange(1.,14./5.,2./5.)*math.pi
    psi.append(inc)
    psi = np.asarray(psi,dtype = np.float64).flatten()
    psi = np.insert(psi,[0,5],[0,0])
    return theta,psi


# In[426]:
#
# Prepare the figure
#
fig = plt.figure()
ax = fig.gca(projection='3d')
plt.ion()
plt.show()
while 1:
    node = int(input("gimme number of node!\n"))
    for artist in plt.gca().lines + plt.gca().collections:

        artist.remove()
    #
    # Find shell numbers
    #
    neigh = dset[node][:,:,0]
    neigh = np.append(neigh,[node])
    pairNum = np.floor(neigh/16.)
    shellNum = 2*pairNum + ((neigh-pairNum*16) > 5.)
    
    nAr = ((shellNum.flatten()))
    pAr = np.unique((pairNum.flatten()))
    g = math.sqrt((1+math.sqrt(5))/2)
    lamb = [math.sqrt(math.sqrt(5)/3),1]
    
    gAr = g**nAr
    gAr = np.unique(gAr)
    
    thetaD,psiD = initializeDodecahedron()
    thetaIc,psiIc = initializeIcosahedron()
    ptsIc = []
    ptsD = []
    pts1= []

    #
    # get each shell's coordinates
    #
    for i,j in enumerate(np.unique(nAr)):
        if np.mod(j,2):
                inc = gAr[i]*lamb[1]*np.array([np.sin(thetaD)*np.cos(psiD),np.sin(thetaD)*np.sin(psiD),np.cos(thetaD)]).T
                ptsD.append(inc)
        else:
                inc = gAr[i]*lamb[0]*np.array([np.sin(thetaIc)*np.cos(psiIc),np.sin(thetaIc)*np.sin(psiIc),np.cos(thetaIc)]).T
                ptsIc.append(inc)
        pts1.append(inc)
    ptsIc = np.asarray(ptsIc,dtype = np.float64)
    ptsD = np.asarray(ptsD,dtype = np.float64)


    #
    # plot interacting triads, triad by triad, for a given vector k[node]
    #
    k_main = k[node,:]
    lims = (-np.amax(gAr*lamb[0]), np.amax(gAr*lamb[0]))
    x = 0.
    y = 0.
    z = 0.
    vX = k_main[0]
    vY = k_main[1]
    vZ = k_main[2]
    conj = np.asarray([1.,-1.])
    ax.quiver(x,y,z,vX,vY,vZ,color = 'r')
    for i in dset[node][:,:,:]:
        k_pair = -conj[i[:,1],None]*k[i[:,0],:]
        x = np.asarray([0, k_pair[0,0]])
        y = np.asarray([0, k_pair[0,1]])
        z = np.asarray([0, k_pair[0,2]])
        vX = k_pair[:,0]
        vY = k_pair[:,1]
        vZ = k_pair[:,2]
        q = ax.quiver(x.T,y.T,z.T,vX.T,vY.T,vZ.T,color = 'k')
        q.set_headlength = 0.1
        
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_zlim(lims)
    for i in ["x", "y", "z"]:
        eval("ax.set_{:s}label('{:s}')".format(i, i))
    for i,j in enumerate(gAr):
    	pts = pts1[i]
    	hull = ConvexHull(pts,qhull_options = 'QJ')
    	for s in range(0,hull.simplices.shape[0]):
    		#s = np.append(s, s[0])  # Here we cycle back to the first coordinate
    		vtx = np.asarray([pts.T[0][hull.simplices[s]], pts.T[1][hull.simplices[s]], pts.T[2][hull.simplices[s]]]).T
    		tri = a3.art3d.Poly3DCollection([vtx])
    		tri.set_edgecolor([93/256,167/256,228/256])
    		tri.set_alpha(0.1)
    		ax.add_collection3d(tri)

    fig.canvas.draw()





