# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 12:56:52 2019

@author: hindesa
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

def EmbedAlg(paramList, guess, F):
    X = []
    for par in paramList:
        try:
            # Optimize F(X,par) = 0
            xval = newton(F, guess, fprime=None, args=(par,), tol=1E-7, maxiter=100)
            # Record solution
            X.append(xval)
            guess = xval
        except RuntimeError:
            #If newton fails, return truncated list of pars with corresponding xvals
            return paramList[:len(X)], X
    # return list of pars and xvals
    return paramList, X

def F(x, lab):
    return lab*x - x**3

labs = np.linspace(5, -5, 200)
C1, X1 = EmbedAlg(labs, 2, F)
C2, X2 = EmbedAlg(labs, -2, F)
C3, X3 = EmbedAlg(labs, 0, F)

#plt.plot(C1, X1)
#plt.plot(C2, X2)
#plt.plot(C3, X3)
k = 8
def budworm(x, k):
    dxdt = r*x*(1-x/k) - x**2/(1+x**2)
    return dxdt

def budLeft(x,r,k):
    y = r*(1-x/k)
    return y

def budRight(x):
    y = x/(1+x**2)
    return y

r = 0.56
k = 8
xs = np.linspace(0,10,200)
ys1 = budLeft(xs,r,k)
ys2 = budRight(xs)
#plt.plot(xs,ys1)
#plt.plot(xs,ys2)

ks = np.linspace(2,10,200)
ks2 = np.linspace(14,6,200)

ks3 = np.linspace(6.1,9.9, 200)

K1, X1 = EmbedAlg(ks, 0.5, budworm)
K2, X2 = EmbedAlg(ks2, 10, budworm)
K3, X3 = EmbedAlg(ks3, 10, budworm)
plt.plot(K1,X1)
plt.plot(K2,X2)
#plt.plot(K3,X3, ls='--')

def px(x):
    return 1/(x*(1+x**2))

def betaLine(px):
    a = 1-px
    b = 2-px
    return a/b

def dSurf(px,beta):
    val = (1-px-beta*px)/beta
    return val

xs = np.linspace(3,11)
pxs = np.linspace(0.1,1.5,1000)
bs = betaLine(pxs)
p = plt.figure()
plt.plot(pxs, bs)

PX1 = px(np.array(X1))
BX1 = betaLine(PX1)

plt.plot(PX1,BX1)

PX2 = px(np.array(X2))
BX2 = betaLine(PX2)

plt.plot(PX2,BX2)

p = plt.figure()
pxs = np.linspace(0.1,1.5,1000)
beta = np.linspace(0,1,1000)
PX, B = np.meshgrid(pxs,beta)
DX = dSurf(PX,B)
