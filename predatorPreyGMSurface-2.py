# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# X = prey, Y = predators

# Samples
n = 100

#S(X) prey births, G(X,Y) predation
#M(Y) predator mortality

s_x =  np.linspace(0,0.5,n)
m_y = np.linspace(0.01,2,n)
g_y = 1



def saddleCondGX(x,y):
    #s_x = x
    #m_y = y
    return x - (x*g_y)/y


X,Y = np.meshgrid(s_x, m_y)

zs = np.array(saddleCondGX(np.ravel(X), np.ravel(Y)))
Z = zs.reshape(X.shape)

fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, Z)
ax2.set_xlabel('s_x')
ax2.set_ylabel('m_y')
ax2.set_zlabel('g_x')
plt.show()




