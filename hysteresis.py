# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 13:43:31 2019

@author: Adrian Hindes
Following http://systems-sciences.uni-graz.at/etextbook/sw2/crittrans.html
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sympy as sp

fig = plt.figure()

# define model parameters
h = 1.; b = 1.; p = 12; r=0.6

def hill(x,p):
    y = x**p/(x**p+h**p)
    return y


def model(x,t,a):
    y = a-b*x+r*hill(x,p)
    return y

def deriv(x):
    dgdx = -b+(r*p*x**(p-1)*h**p)/((x**p + h**p)**2)
    return dgdx


minA = 0.4
maxA = 0.9

aS = np.arange(minA,maxA,0.1)

XS = []
x0 = 2
t = np.linspace(0,10,num=10000)
for a in aS:
    x = integrate.odeint(model,x0,t,args=(a,))
    plt.plot(t,x,label='a = '+str(round(10*a)/10))
    XS.append(x)
    
plt.xlabel('t')
plt.ylabel('x')
plt.grid()
plt.legend(loc='best')

fig2 = plt.figure()
n = 1000
A = np.linspace(minA,maxA,n)
XA = np.zeros(n)

for j in range(n):
    x = x0;
    xs = integrate.odeint(model,x0,t,args=(A[j],))
    xeq = xs[-1][0]
    XA[j] = xeq
    
plt.plot(A,XA)