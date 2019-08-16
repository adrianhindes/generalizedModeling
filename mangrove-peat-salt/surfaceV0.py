# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:22:22 2019

@author: hindesa
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# M = mangrove biomass (kg), P = peat soil elevation (mm)
# S = soil salinity (ppm?)

# Linspace range
n = 1

# ---------------
# Timescales
# ---------------
# 1/years turnover rates
# total guesses at the moment
alphaM = 1.
alphaP = 1.
alphaS = 1.



# ----------------
# Beta parameters
# ----------------
# Must add to 1
# Mangrove gain
betaP = random.uniform(0,0.5,n) # from propagules & established saplings
betaG = 1-betaP # from endogenous growth of existing trees

# Mangrove loss
betaD = 0.3
betaS = 0.5
betaL = 1-betaD-betaS

# Peat gain
betaA = 0.6
betaR = 0.1
betaV = 1 - betaA - betaR

#Peat loss
betaE = 0.5
betaSB = 1 - betaE

# ----------------------
# Elasticity parameters
# ----------------------
#hydP = random.uniform(-2.0, 0.0, n)
# Mangroves
propM = 1 
propS = -1
growM = 1

drownHyd = 2
drownM = 1 

stressM = 1
stressS = 2

littM = 1.5

# Peat soils
accSed = 1 
sedHyd = 2
accM = 1

retLitt = 1
retHyd = 1

volGrow = 1
volP = 1

eroM = 1

subsM = 1
subsHyd = 1
subsP = 0.5

# Salinity
inM = 1.0 
inS = 0.5

outS = 1

# Remap parameters for surface plot
k = 1000
betaD = np.linspace(0,1, k)
betaL = 0.1 # minimal loss from leaf litter
betaE = np.linspace(0,1, k)


def saddleHyd(betaD,betaE):
    betaS = 1-betaD-betaL
    betaSB = 1-betaE
    
    hyd = -((betaP*betaV*outS-betaP*betaV*inS)*volP+(betaP*betaSB*inS-betaP*betaSB*outS)*\
         subsP-betaS*inM*stressS+betaP*inM*propS)/\
        ((betaD*betaV*drownHyd*growM*outS\
        -betaD*betaV*drownHyd*growM*inS)*volGrow+(betaD*betaSB*drownHyd*inS\
        -betaD*betaSB*drownHyd*outS)*subsM+(betaP*betaSB*inS-betaP*betaSB*outS)\
        *subsHyd+(accSed*betaA*betaP*outS-accSed*betaA*betaP*inS)*sedHyd\
        +(betaD*betaR*drownHyd*littM*outS-betaD*betaR*drownHyd*inS*littM)\
        *retLitt+(betaP*betaR*outS-betaP*betaR*inS)*retHyd+(accM*betaA*betaD*drownHyd-betaD*betaE*drownHyd*eroM)\
        *outS+(betaD*betaE*drownHyd*eroM-accM*betaA*betaD*drownHyd)*inS)
    return hyd


X,Y = np.meshgrid(betaD,betaE)

zs = np.array(saddleHyd(np.ravel(X), np.ravel(Y)))
Z = zs.reshape(X.shape)

fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

ax2.plot_surface(Y, X, Z)
ax2.set_xlabel('betaD')
ax2.set_ylabel('betaE')
ax2.set_zlabel('hydP')
plt.show()

