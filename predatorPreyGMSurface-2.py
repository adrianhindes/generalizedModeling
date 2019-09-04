# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('text', usetex=True)

# X = prey, Y = predators

# Samples
n = 100

#S(X) prey births, F(X,Y) predation
#D(X) Prey mortality, M(Y) predator mortality


m_y = 1
d_x = 1
f_x = 2

beta = np.linspace(0,1,n)
f_y = np.linspace(0,0.8,n)


def saddleCondGX(beta,f_y):
    beta2 = 1-beta
    
    return beta + 2*beta2*(1-f_y/(f_y-1))

def saddle(beta,f_y):
    beta2 = 1-beta
    a = beta*d_x + beta2*f_x
    b = f_y-m_y
    c = beta2*f_x*f_y
    
    res = (a*b-c)/b
    return res



X,Y = np.meshgrid(beta, f_y)

Z = saddleCondGX(X,Y)

fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, Z)
ax2.set_xlabel(r'$\beta$')
ax2.set_ylabel(r'$f_y$')
ax2.set_zlabel(r'$s_x$')
plt.show()




