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
beta_m = 0.2
beta_m2 = 1 - beta_m

# Elasticity parameters
g_m = np.linspace(0, 2, n) # elasticity of growth w resp mangroves
g_w = 0.2 # elasticity of growth w resp water depth

d_m = 1 # elasticity of endogenous mortality

i_w = 2 # elasticity of inundation induced dieback

#a_m = 1
s_p = 1 # elasticity of subsidence w resp. soil elevation
s_w = 1 # elasticity of subsidence w resp. water depth

q_p = np.linspace(0, 1, n) # elasticity of water displacement w resp soil elevation

# Function for saddle-node bifurcation surface
# w respect to  q_p and g_m

def saddle(gm, qp, iw):
    a = gm-beta_m*d_m
    b = s_p - qp*s_w
    c = qp*(beta_m2*iw-g_w)
    
    a_m = a*b/c

    return a_m

iw0=0.2

X,Y = np.meshgrid(g_m, q_p)

zs = np.array(saddle(np.ravel(X), np.ravel(Y), iw0))
Z = zs.reshape(X.shape)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

l = ax.plot_surface(X,Y,Z, rstride=3, cstride=3)         
ax.set_xlabel('g_m')
ax.set_ylabel('q_p')
ax.set_zlabel('a_m')

axSlider = plt.axes([0.1, 0.05, 0.65, 0.03])
paramSlider = Slider(axSlider, 'i_w', 0.0, 2, valinit=iw0)

def update(val): 
    iw = paramSlider.val 
    
    ax.clear()
    Z = saddle(X, Y, iw)
    
    l=ax.plot_surface(X,Y,Z,rstride=3, cstride=3)
    
    ax.set_xlabel('g_m')
    ax.set_ylabel('q_p')
    ax.set_zlabel('a_m')

    
    fig.canvas.draw_idle() 

paramSlider.on_changed(update)

plt.show()


