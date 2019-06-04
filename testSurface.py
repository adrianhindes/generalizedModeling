# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:03:45 2019

@author: hindesa
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

n=100
def f(x,y):
    e1 = x-2*y+1
    e2 = 1-y
    z = e1/e2
    return z

x = np.linspace(0,1,n)
y = np.linspace(0,1,n)
X,Y = np.meshgrid(x,y)

Z = f(X,Y)
fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, Z)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')
plt.show()