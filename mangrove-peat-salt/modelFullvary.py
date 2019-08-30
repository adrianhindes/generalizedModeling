# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
import numpy as np
from numpy import random
from numpy import linalg as LA
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import stats
import pandas as pd
import seaborn as sns
from math import isnan
from jacobian import computeJac

# M = mangrove biomass (kg), P = peat soil elevation (mm)
# S = soil salinity (ppm?)

# Linspace range
n = 10000

# ---------------
# Timescales
# ---------------
# 1/years turnover rates
# total guesses at the moment
alphaM = random.uniform(1/6,1/6,n)
alphaP = random.uniform(1/3,1/3,n)
alphaS = random.uniform(3,3,n)

paramLabels = {'alphaM','alphaP','alphaS',
               'betaP','betaD','betaL',
        'betaA','betaR','betaE',
        'hydP','propM','propS','growM',
        'drownHyd','drownM','stressM',
        'stressS','littM','accSed',
        'sedHyd','accM','retLitt','retHyd',
        'volGrow','volP','eroM','subsM',
        'subsHyd','subsP','inS','inM','outS'}

# ----------------
# Beta parameters
# ----------------
# Must add to 1
# Mangrove gain
betaP = random.uniform(0,1,n) # from propagules & established saplings
betaG = 1-betaP # from endogenous growth of existing trees

# Mangrove (biomass) loss
betaL = random.uniform(0,1,n) # litter fall
betaD = random.uniform(0,1-betaL, n) # drowning
betaS = 1-betaL-betaD # salt stress


# Peat gain
betaR = random.uniform(0,1,n) # retained litterfall
betaA = random.uniform(0, 1-betaR, n) #sediment accretion
betaV = 1 - betaA - betaR #volume incrcease

#Peat loss
betaE = random.uniform(0,1,n)
betaSB = 1 - betaE

# ----------------------
# Elasticity parameters
# ----------------------

# Default range for unknown but likely linear values
r0 = 0.5
r1 = 2


hydP = random.uniform(-2.0, 0, n)
# Mangroves
propM = random.uniform(r0, r1, n) 
propS = random.uniform(-1, 0.0, n)
growM = random.uniform(r0, r1, n)

drownHyd = random.uniform(0.0, 5.0, n)
drownM = random.uniform(r0, r1, n) 

stressM = random.uniform(r0, r1, n)
stressS = random.uniform(0.0, 2.0, n)

littM = random.uniform(r0, r1, n)

# Peat soils
accSed = random.uniform(r0, r1, n) 
sedHyd = random.uniform(0.5, 4, n)
accM = random.uniform(r0, r1, n)

retLitt = random.uniform(r0, r1, n)
retHyd = random.uniform(-2.0, 0.0, n)

volGrow = random.uniform(r0, r1, n)
volP = random.uniform(r0, r1, n)

eroM = random.uniform(-3,0, n)

subsM = random.uniform(r0, r1, n)
subsHyd = random.uniform(r0, r1, n)
subsP = random.uniform(r0, r1, n)

# Salinity
inM = random.uniform(r0,r1,n)
inS = random.uniform(r0, 1, n)

outS = random.uniform(r0, 1,n)

def stability(eigs):
    # take in vector of eigenvalues
    # tell if system stable or not (stable if all eigenvalues are less than zero)
    reals = np.real(eigs)
    if max(reals) < 0:
        result = 1
    else:
        result = 0
    return result
# Construct dataframe to track parameters and associated eigenvalues
# Parameters that are varying
data = {'alphaM':alphaM,'alphaP':alphaP,'alphaS':alphaS,
        'betaP':betaP,'betaD':betaD,'betaL':betaL,
        'betaA':betaA,'betaR':betaR,'betaE':betaE,
        'hydP':hydP,'propM':propM,'propS':propS,'growM':growM,
        'drownHyd':drownHyd,'drownM':drownM,'stressM':stressM,
        'stressS':stressS,'littM':littM,'accSed':accSed,
        'sedHyd':sedHyd,'accM':accM,'retLitt':retLitt,'retHyd':retHyd,
        'volGrow':volGrow,'volP':volP,'eroM':eroM,'subsM':subsM,
        'subsHyd':subsHyd,'subsP':subsP,'inS':inS,'inM':inM,'outS':outS}

eigs = [] #eigenvalue triplets
eigsV = [] #eigenvector triplets
stab = [] #stability (0,1)
determ = [] # determinant of Jacobian
for j in tqdm(range(n)):
    
    dataJ = {k:v[j] for (k,v) in data.items()}
    jac = computeJac(dataJ)
   
    w, v = LA.eig(jac)
    det = LA.det(jac)
    
    eigs.append(np.real(np.max(w)))
    eigsV.append(v)
    stab.append(stability(w))
    determ.append(det)

# Compute correlations
corrs = {k:(np.corrcoef(data[k],stab)[0,1]) for (k,v) in data.items() }

#Sort out only the big correlations
numPlot = 15 # number of variables to plot
remNum = len(corrs.items()) - numPlot #number of variables to remove
absCorrs = {k: np.abs(v) for (k,v) in corrs.items() if not isnan(v)}
corrsSorted = sorted(absCorrs.items(), key=lambda x: x[1], reverse=True)
delN = len(corrsSorted) - numPlot
del corrsSorted[-delN:]

bigCorrs = {k: corrs[k] for (k,v) in corrsSorted}

p2 = plt.bar(range(len(bigCorrs)), bigCorrs.values())

plt.ylabel('Correlation Coefficient')
plt.xlabel('Parameter')
plt.xticks(range(len(bigCorrs)), list(bigCorrs.keys()), rotation=70)
plt.show()

posStable = set() # parameters for which higher values => greater stability
negStable = set() # parameters for which higher values => lower stability

for (param,corr) in bigCorrs.items():
    if corr > 0 :
        posStable.add(param)
    else:
        negStable.add(param)
        

# For max stability
maxStable = {}
minStable = {}

for param in posStable:
    maxStable[param] = max(data[param])
    minStable[param] = min(data[param])
    
for param in negStable:
    maxStable[param] = min(data[param])
    minStable[param] = max(data[param])
    
