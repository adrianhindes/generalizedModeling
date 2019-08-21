# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:39:53 2019

@author: Adrian
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 100

def cubic(x,a):
    y = x**3+a
    return y

x = np.linspace(-10,10,n)
a = np.linspace(0,1,n)

c = 1

plt.plot(x,cubic(x,c))

X,A = np.meshgrid(x,a)
Y = cubic(X,A)

fig1=plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.plot_surface(X, A, Y)
ax1.set_xlabel('x')
ax1.set_ylabel('a')
ax1.set_zlabel('y')
plt.show()