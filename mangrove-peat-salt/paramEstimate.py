# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:52:10 2019

@author: Adrian
Parameter estimations
"""
import numpy as np
#Lovelock 2017 propS 
A = 22149609
sigA = 22954957
c = -0.2035
sigC = 0.0208

def elas(salt):
    # compute elasticity of regression given an assumed soil porewater salinity level
    # corresponding to a fixed point    
    y = A*np.exp(c*salt)
    dyds = c*A*np.exp(c*salt)
    elas = (y/salt)*dyds
    return elas
def elasPlus(salt):
    y = (A+sigA)*np.exp((c+sigC)*salt)
    dyds = (c+sigC)*(A+sigA)*np.exp((c+sigC)*salt)
    elas = (y/salt)*dyds
    return elas
def elasMin(salt):
    y = (A-sigA)*np.exp((c-sigC)*salt)
    dyds = (c-sigC)*(A-sigA)*np.exp((c-sigC)*salt)
    elas = (y/salt)*dyds
    return elas

