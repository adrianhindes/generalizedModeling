# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np

# M = mangrove biomass(?), P = peat soil elevation

# Linspace range
n = 100
# Timescales


# Beta parameters
# beta_m = proportion of mangrove loss from natural mortality
# beta_m2 = proportion of mangrove loss from coastal inundation
branch = dict()
branch['betam'] = 0.2

# Elasticity parameters
elast = dict()
elast['g_m'] = 1
elast['g_w'] = 0.2
elast['d_m'] = 1 # elasticity of endogenous mortality

elast['i_w'] = 2 # elasticity of inundation induced dieback

elast['a_m'] = 1
elast['s_p'] = 1 # elasticity of subsidence w resp. soil elevation
elast['s_w'] = 1 # elasticity of subsidence w resp. water depth

elast['q_p'] = 0.5 # elasticity of water displacement w resp soil elevation

# Function for saddle-node bifurcation surface
# w respect to  q_p and g_m

def saddle(branch, elast):
    betam = branch['betam']
    betam2 = 1-betam
    # Mangrove gm parameters
    g_m = elast['g_m']
    g_w = elast['g_w']
    d_m = elast['d_m']
    i_w = elast['i_w']
    
    # Peat soil gm parameters
    a_m = elast['a_m']
    s_p = elast['s_p']
    s_w = elast['s_w']
    
    # Water gm parameter
    q_p = elast['q_p']
    
    a = g_m-betam*d_m
    b = s_p-q_p*s_w
    c = a_m*q_p*g_w
    d = a_m*q_p*betam2
    
    iw = (a*b+c)/d
    return iw

beta_m0=0.2
x = np.linspace(0,2,n)
y = np.linspace(0,1,n)

X,Y = np.meshgrid(x,y)

elast['g_m'] = x
elast['q_p'] = y
branch['beta_m0'] = 0.2

Z = np.array(saddle(branch, elast))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

l = ax.plot_surface(X,Y,Z, rstride=3, cstride=3)         
ax.set_xlabel('g_m')
ax.set_ylabel('q_p')
ax.set_zlabel('i_w')


axBeta = plt.axes([0.1, 0.05, 0.65, 0.02])
sBeta = Slider(axBeta, 'beta_m', 0.0, 1.0, valinit=beta_m0)

def update(val): 
    beta_m = sBeta.val 
    ax.clear()
    Z = saddle(X, Y, beta_m)
    
    l=ax.plot_surface(X,Y,Z,rstride=3, cstride=3)
    
    ax.set_xlabel('g_m')
    ax.set_ylabel('q_p')
    ax.set_zlabel('i_w')
    

    
    fig.canvas.draw_idle() 

sBeta.on_changed(update)

plt.show()


