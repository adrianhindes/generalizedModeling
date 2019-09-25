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
            xval = newton(F, guess, fprime=Fprime, args=(par,), tol=1E-8, maxiter=100)
            # Record solution
            X.append(xval)
            guess = xval
        except RuntimeError:
            #If newton fails, return truncated list of pars with corresponding xvals
            return paramList[:len(X)], X
    # return list of pars and xvals
    return paramList, X

def F(x, lab):
    return lab + x - x**3

def Fprime(x, lab):
    return 1-3*x**2

labs = np.linspace(3, -1/np.sqrt(3)+0.02, 200)
labs2 = np.linspace(-3,1/np.sqrt(3)-0.01, 200)
labs3 = np.linspace(-1/np.sqrt(3)+0.2,1/np.sqrt(3)-0.2,200)
C1, X1 = EmbedAlg(labs, 1, F)
C2, X2 = EmbedAlg(labs2, -1, F)
C3, X3 = EmbedAlg(labs3, 0, F)

plt.plot(C1, X1)
plt.plot(C2, X2)
plt.plot(C3, X3, ls='--')
