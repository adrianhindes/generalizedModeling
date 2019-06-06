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





# Beta parameters
# beta_m = proportion of mangrove loss from natural mortality
# beta_m2 = proportion of mangrove loss from coastal inundation
beta_m = np.linspace(0,1, n)
beta_m2 = 1 - beta_m

# Elasticity parameters
g_m = 1.5 # elasticity of growth w resp mangroves
g_w = 0.5 # elasticity of growth w resp water depth

d_m = 1 # elasticity of endogenous mortality

#i_w = 2 # elasticity of inundation induced dieback

a_m = 2 # elasticity of peat accumulation w resp. mangroves
s_p = 1 # elasticity of subsidence w resp. soil elevation
s_w = np.linspace(0,2,n) # elasticity of subsidence w resp. water depth

q_p = 0.5 # elasticity of water displacement w resp soil elevation

# Function for saddle-node bifurcation surface
# w respect to  q_p and g_m

def saddle(beta_m2, s_w, g_m):
    beta_m = 1-beta_m2
    
    a = g_m-beta_m*d_m
    b = s_p-q_p*s_w
    c = a_m*q_p*g_w
    d = a_m*q_p*beta_m2
    
    iw = (a*b+c)/d
    return iw

g_m0=1

X,Y = np.meshgrid(beta_m2, s_w)

Z = saddle(X, Y, g_m)
Z[Z >= 10]=10
Z[Z <= -10]=-10

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

l = ax.plot_surface(X,Y,Z, rstride=3, cstride=3)         
ax.set_xlabel('Prop. inundation deaths')
ax.set_ylabel('Peat subsidence elasticity (water)')
ax.set_zlabel('Inundation elasticity')

axSlider = plt.axes([0.1, 0.05, 0.65, 0.03])
growthSlider = Slider(axSlider, 'Growth', 0.0, 3.0, valinit=g_m0)

def update(val): 
    g_m = growthSlider.val 
    ax.clear()
    Z = saddle(X, Y, g_m)
    Z[Z >= 10]=10
    Z[Z <= -10]=-10

    l=ax.plot_surface(X,Y,Z,rstride=3, cstride=3)
    
    ax.set_xlabel('Prop. inundation deaths')
    ax.set_ylabel('Peat subsidence elasticity (water)')
    ax.set_zlabel('Inundation elasticity')

    
    fig.canvas.draw_idle() 

growthSlider.on_changed(update)

plt.show()


