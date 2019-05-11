# -*- coding: utf-8 -*-
"""
Created on Sat May 11 15:15:56 2019

@author: hindesa
Simple bifurcation diagram test
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import integrate

def dX_dt(x, r, t=0):
    """ Return parabolic growth rate of x with param r"""
    return r*x + x**3 - x**5

points = 1000
fig = plt.figure()
ax = fig.gca(projection='3d')

x = np.linspace(-10, 10, points)
r = np.linspace(-10, 10, points)

X,R, = np.meshgrid(x,r)
Z = dX_dt(X,R)

surf = ax.plot_surface(X, R, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()

# Numerically calculate values close to equilibrium
tol = 0.05 #error tolerance
eqIndex = [(i,j) for i in range(points) for j in range(points) if np.hypot(Z[i,j], 0.0) <= tol ]

#Recover values
Xeq = [X[point[0],point[1]] for point in eqIndex]
Req = [R[point[0],point[1]] for point in eqIndex]

fig2 = plt.figure()
plt.scatter(Req,Xeq)
plt.xlabel('Parameter value r')
plt.ylabel('Equilibrium point x')
plt.title('Ostensible Bifurcation plot of dx/dt = x^2 +r')
plt.show()