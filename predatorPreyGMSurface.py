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
# Timescales
# Fix timescales to be equivalent to foxes and rabbits
# avg fox lifetime 4 years
# avg rabbit lifetime 1.5 years
alphax = 1/4
alphay = 1/1.5

# Beta parameters
beta1 = np.random.uniform(low=0.0, high = 1.0, size=n)
beta2 = 1 - beta1

# Elasticity parameters
# Only vary elasticities of prey deaths
f_x = np.random.uniform(low=0.0, high = 2.0, size=n) # Predation elasticity w resp. x
m_x = f_x # Mortality elasticiy of prey
f_y = 2 # quadratic predation rate for predators
m_y = 1 # Mortality elasticiy of predators


def saddleCond(x,y):
    # Smaller calculations
    # beta = x
    # f_x = y
    x2 = 1-x
    a = y*(x+x2)
    b = f_y-m_y
    c = x2*y*f_y
    s_x = (a*b-c)/b
    return s_x

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#X,Y = np.meshgrid(beta1,f_x)
#
#zs = np.array(saddleCond(np.ravel(X), np.ravel(Y)))
#Z = zs.reshape(X.shape)
#
#ax.plot_surface(X, Y, Z)
#ax.set_xlabel('Relative contribution parameter for deaths')
#ax.set_ylabel('Elasticity of predation rate w resp. prey')
#ax.set_zlabel('Elasticiy of prey birth rate')
#plt.show()


# Remap parameters for Hopf surface
alpha_r = np.linspace(0,1, n)
beta = np.linspace(0,1, n)


def hopfCond(alphaR, beta):
    #Hopf condition in terms of alphaR, beta
    s_x = 2
    d_x = 1
    f_y = 1
    
    e1 = alphaR*(f_y-m_y)
    e2 = s_x - beta*d_x
    
    fx = (1/(1-beta))*(e1+e2)

    return fx

X,Y = np.meshgrid(alpha_r, beta)

zs = np.array(hopfCond(np.ravel(X), np.ravel(Y)))
Z = zs.reshape(X.shape)

fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(Y, X, Z)
plt.xlim(1,0)
ax2.set_xlabel('alpha_r')
ax2.set_ylabel('beta')
ax2.set_zlabel('fx')
plt.show()



