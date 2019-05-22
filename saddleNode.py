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
from pynamical import simulate, bifurcation_plot
from numba import jit

@jit(cache=True, nopython=True) # pragma: no cover 
def map2(x, r):
    """
    Define the equation for the cubic map.
    
    Arguments
    ---------
    x: float
        current population value at time t
    rate: float
        growth rate parameter values
    
    Returns
    -------
    float
         result of map at time t+1
    """
    
    return r - x - np.exp(x)


x = simulate(model=map2, num_gens = 100, initial_pop=1, rate_min=0, rate_max=5, num_rates=2000, num_discard=100)
bifurcation_plot(x, xmin=-10, xmax=10, ymin=-10, ymax=10)