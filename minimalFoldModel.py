# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 13:43:31 2019

@author: Excalibur
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sympy as sp

a = 0
b = 0.5
h = 2
p = 2
r = 10



def hill(x):
    p=2
    y = x**p/(x**p+h**p) 
    return y

x = sp.Symbol('x')
rs = np.linspace(0,1,50)
states = []

def stables(r):
    modFun = a-b*x+r*hill(x)
    stableStates = sp.solve(sp.Eq(0,modFun),x)
    return stableStates

for r in rs:
    states.append(stables(r))    


def minMod(t,x,r):
    dxdt = a-b*x+r*hill(x)
    return dxdt

def ode(t, x):
    r = t/1000
    dxdt = minMod(t,x,r)
    return dxdt

x0 = [5]
n = 100
tEnd = 1000
t = np.linspace(0,tEnd,n)


sol = integrate.solve_ivp(fun=ode, t_span=[0, tEnd], y0=x0, t_eval=t)
xs = sol.y[0]

plt.plot(t, xs)