# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:25:24 2019

@author: hindesa
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


k = 100
r = 0.5

def logistic(N,t):
    return r*N*(1-(1/k)*N)
    
ts = np.linspace(0,20,1000)
n0 = 1
ys = odeint(logistic, n0, ts)
ys = np.array(ys).flatten()

plt.xlabel('Time')
plt.ylabel('Population (N)')
plt.title('Logistic Population Growth, r = 0.5, K=100')
plt.plot(ts,ys)