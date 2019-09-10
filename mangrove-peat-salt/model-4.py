                    # -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
from numpy import random
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np
from tqdm import tqdm
# M = mangrove biomass (kg), P = peat soil elevation (mm)
# S = soil salinity (ppm?)

# Linspace range
n = 10000

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
hydP = random.uniform(-2.0, 0.0, n)
# Mangroves
propM = 1 
propS = -5
growM = 1

drownHyd = random.uniform(0.0, 5.0, n)
drownM = 1 

stressM = 1
stressS = random.uniform(0.0, 2.0, n)

littM = random.uniform(1.0, 2.0, n)

# Peat soils
accSed = 1 
sedHyd = random.uniform(0.0, 2.0, n)
accM = 1

retLitt = 1
retHyd = random.uniform(-2.0, 0.0, n)

volGrow = 1
volP = 1

eroM = 1

subsM = 1
subsHyd = random.uniform(0.0, 1.0, n)
subsP = random.uniform(0.5, 1.5, n)

# Salinity
inM = 1.0 
inS = random.uniform(0.0, 1.0, n)

outS = 1

# Define Jacobian matrix elements
# Note syntax here is not strictly correct - dmdm == (dm/dt)/dm

hydP = np.linspace(-2.0,0.0,n)


dmdm = betaP*propM +betaG*growM -betaS*stressM -betaD*drownM -betaL*littM
dmdp = -betaD*hydP*drownHyd
dmds = betaP*propS -betaS*stressS

dpdm = betaA*accM +betaR*retLitt*littM + betaV*volGrow*growM\
        -betaE*eroM -betaSB*subsM
        
dpdp = hydP*(betaA*accSed*sedHyd + betaR*retHyd-betaSB*subsHyd)\
        +betaV*volP-betaSB*subsP
        
dpds = np.zeros(n)

dsdm = [inM]*n
dsdp = np.zeros(n)
dsds = inS - outS


# Get eigenvalues
eigens = []
R1 = [dmdm, dmdp, dmds]
R2 = [dpdm, dpdp, dpds]
R3 = [dsdm, dsdp, dsds]


for j in tqdm(range(n)):
    jac =  np.array([ [dmdm[j], dmdp[j], dmds[j]],
                      [dpdm[j], dpdp[j], dpds[j]],
                      [dsdm[j], dsdp[j], dsds[j]]
                    ])
    w, v = LA.eig(jac)
    eigens.append(w)
    
    
def stability(eigs):
    # take in vector of eigenvalues
    # tell if system stable or not (stable if all eigenvalues are less than zero)
    reals = np.real(eigs)
    if max(reals) < 0:
        result = 1
    else:
        result = 0
    return result

stabs = list(map(stability, eigens))
numStable = 0
for j in range(len(stabs)):
    numStable += stabs[j]


maxReals = np.array([np.real(max(v)) for v in eigens])


plt.scatter(hydP, maxReals)