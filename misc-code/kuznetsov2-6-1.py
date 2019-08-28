# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:47:34 2019

@author: hindesa

Lotka-Volterra specific case
Following SciPy Tutorial:
    https://scipy-cookbook.readthedocs.io/items/LoktaVolterraTutorial.html
"""

import numpy as np
import pylab as plt

def dX_dt(A, t=0):
    """ Return the growth rate of prey and predator population """
    x = A[0]
    y = A[1]
    matrix = np.array([x + 2*y,
                    -x-y])
    return matrix

n = 25
a = -5.0
b = 5.0
x = np.linspace(a, b, n)
y = np.linspace(a, b, n)

X1,Y1 = np.meshgrid(x,y)

t = 0

u, v = np.zeros(X1.shape), np.zeros(Y1.shape)

lenI, lenJ = Y1.shape
# Calculate quiver values
for i in range(lenI):
    for j in range(lenJ):
        x = X1[i,j]
        y = Y1[i,j]
        deriv = dX_dt([x,y], t)
        u[i,j] = deriv[0]
        v[i,j] = deriv[0]

f1 = plt.figure()
Q = plt.quiver(X1,Y1, u, v, color='r')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.xlim([a, b])
plt.ylim([a, b])
plt.show()
