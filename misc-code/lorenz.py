# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 16:31:32 2019

@author: Excalibur
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Lorenz parameters and init conditions
sigma, beta, rho = 10, 2.667, 28
u0, v0, w0 = 0, 1, 1.05

#Time for simulation and steps
tmax, n = 100, 10000

def lorenz(X, t, sigma, beta, rho):
    u,v,w = X
    up = -sigma*(u-v)
    vp = rho*u - v - u*w
    wp = -beta*w+u*v
    return up, vp, wp

# Integrate
t = np.linspace(0,tmax,n)
f = odeint(lorenz, (u0,v0,w0),t,args=(sigma,beta,rho))
x, y, z = f.T

fig = plt.figure()
ax = fig.gca(projection='3d')
s=10
c = np.linspace(0,1,n)
for i in range(0,n-s,s):
    ax.plot(x[i:i+s+1], y[i:i+s+1], z[i:i+s+1], color=(1,c[i],0),alpha=0.4)

f2 = odeint(lorenz, (u0+0.1,v0,w0),t,args=(sigma,beta,rho))

for i in range(0,n-s,s):
    ax.plot(x[i:i+s+1], y[i:i+s+1], z[i:i+s+1], color=(c[i],0,0),alpha=0.4)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$z$')

plt.show()