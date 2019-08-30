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

from parameterDefaults import *
startParams = [betas, mangs, peats, salts]
from droughtParams import *
endParams = [betas, mangs, peats, salts]
# M = mangrove biomass (kg), P = peat soil elevation (mm)
# S = soil salinity (ppm?)

# Linspace range
n = 10000



# ---------------
# Timescales
# ---------------
# 1/years turnover rates
# total guesses at the moment
alphaM = 1/6.
alphaP = 1/3.
alphaS = 6.


# ----------------
# Beta parameters
# ----------------
betas = {}
for key in set(startParams[0]):
    betas[key] = (startParams[0][key],endParams[0][key])
    
# Must add to 1
# Mangrove gain
betaP = np.linspace(betas['betaP'][0], betas['betaP'][1], n) # from propagules & established saplings
betaG = 1-betaP # from endogenous growth of existing trees

# Mangrove (biomass) loss
betaL = np.linspace(betas['betaL'][0], betas['betaL'][1], n) # litter fall
betaD = np.linspace(betas['betaD'][0], betas['betaD'][1], n)
betaS = 1-betaL-betaD # salt stress


# Peat gain
betaR = np.linspace(betas['betaR'][0], betas['betaR'][1], n) # retained litterfall
betaA = np.linspace(betas['betaA'][0], betas['betaA'][1], n) #sediment accretion
betaV = 1 - betaA - betaR #volume incrcease

#Peat loss
betaE = np.linspace(betas['betaE'][0], betas['betaE'][1], n)
betaSB = 1 - betaE

# ----------------------
# Elasticity parameters
# ----------------------
# Mangroves

mangs = {}
for key in set(startParams[1]):
    mangs[key] = (startParams[1][key],endParams[1][key])
    
propM = np.linspace(mangs['propM'][0], mangs['propM'][1], n) 
propS = np.linspace(mangs['propS'][0], mangs['propS'][1], n) 
growM = np.linspace(mangs['growM'][0], mangs['growM'][1], n) 

drownHyd = np.linspace(mangs['drownHyd'][0], mangs['drownHyd'][1], n) 
drownM = np.linspace(mangs['drownM'][0], mangs['drownM'][1], n) 

stressM = np.linspace(mangs['stressM'][0], mangs['stressM'][1], n) 
stressS = np.linspace(mangs['stressS'][0], mangs['stressS'][1], n) 

littM = np.linspace(mangs['littM'][0], mangs['littM'][1], n) 

# Peat soils
peats = {}
for key in set(startParams[2]):
    peats[key] = (startParams[2][key],endParams[2][key])
    
accSed = np.linspace(peats['accSed'][0], peats['accSed'][1], n) 
sedHyd = np.linspace(peats['sedHyd'][0], peats['sedHyd'][1], n) 
accM = np.linspace(peats['accM'][0], peats['accM'][1], n) 

retLitt = np.linspace(peats['retLitt'][0], peats['retLitt'][1], n) 
retHyd = np.linspace(peats['retHyd'][0], peats['retHyd'][1], n) 

volGrow = np.linspace(peats['volGrow'][0], peats['volGrow'][1], n) 
volP = np.linspace(peats['volP'][0], peats['volP'][1], n) 

eroM = np.linspace(peats['eroM'][0], peats['eroM'][1], n) 

subsM = np.linspace(peats['subsM'][0], peats['subsM'][1], n) 
subsHyd = np.linspace(peats['subsHyd'][0], peats['subsHyd'][1], n)
subsP = np.linspace(peats['subsP'][0], peats['subsP'][1], n) 

hydP = np.linspace(peats['hydP'][0], peats['hydP'][1], n) 

# Salinity
salts = {}
for key in set(startParams[3]):
    salts[key] = (startParams[3][key],endParams[3][key])
    
    
inM = np.linspace(salts['inM'][0],salts['inM'][1], n)
inS = np.linspace(salts['inS'][0],salts['inS'][1], n)

outS = np.linspace(salts['outS'][0],salts['outS'][1], n)

def stability(eig):
    # take in vector of eigenvalues
    # tell if system stable or not (stable if all eigenvalues are less than zero)
    if np.real(eig) < 0:
        result = 1
    else:
        result = 0
    return result

# Construct dataframe to track parameters and associated eigenvalues
# Parameters that are varying
data = {'betaP':betaP,'betaD':betaD,'betaL':betaL,
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
    
    # Remap parameters
    betaP = data['betaP'][j]
    betaD = data['betaD'][j]
    betaL = data['betaL'][j]
    
    betaA = data['betaA'][j]
    betaR = data['betaR'][j]
    betaE = data['betaE'][j]
    hydP = data['hydP'][j]
    propM = data['propM'][j]
    propS = data['propS'][j]
    growM = data['growM'][j]
    
    drownHyd = data['drownHyd'][j]
    drownM = data['drownM'][j]
    stressM = data['stressM'][j]
    stressS = data['stressS'][j]
    
    littM = data['littM'][j]
    accSed = data['accSed'][j]
    sedHyd = data['sedHyd'][j]
    accM = data['accM'][j]
    retLitt = data['retLitt'][j]
    retHyd = data['retHyd'][j]
    volGrow = data['volGrow'][j]
    volP = data['volP'][j]
    
    eroM = data['eroM'][j]
    subsM = data['subsM'][j]
    subsHyd = data['subsHyd'][j]
    subsP = data['subsP'][j]
    inS = data['inS'][j]
    inM = data['inM'][j]
    outS = data['outS'][j]
    
    # Define Jacobian matrix elements
    # Note syntax here is not strictly correct - dmdm == (dm/dt)/dm
    
    betaG = 1-betaP
    betaS = 1-betaD-betaL
    betaSB = 1-betaE
    betaV = 1-betaA-betaR
    
    dmdm = betaP*propM +betaG*growM -betaS*stressM -betaD*drownM -betaL*littM
    dmdp = -1*betaD*hydP*drownHyd
    dmds = betaP*propS -betaS*stressS

    dpdm = betaA*accM +betaR*retLitt*littM + betaV*volGrow*growM\
            -betaE*eroM -betaSB*subsM
        
    dpdp = hydP*(betaA*accSed*sedHyd + betaR*retHyd-betaSB*subsHyd)\
            +betaV*volP-betaSB*subsP
    dpds = 0
    
    dsdm = inM
    dsdp = 0
    dsds = inS - outS
    
    # alpha paramater array
    alphas = np.array([ [alphaM, 0, 0], [0, alphaP, 0], [0, 0, alphaS]])

    R1 = [dmdm, dmdp, dmds]
    R2 = [dpdm, dpdp, dpds]
    R3 = [dsdm, dsdp, dsds]

    jac0 =  np.array([[dmdm, dmdp, dmds],
                      [dpdm, dpdp, dpds],
                      [dsdm, dsdp, dsds]])
    jac = np.matmul(alphas, jac0)
    w, v = LA.eig(jac)
    det = LA.det(jac)
    
    eigMax = np.max(w)
    eigVMax = v[: , w.argmax()]
    
    eigs.append(eigMax)
    eigsV.append(eigVMax)
    stable = np.real(eigMax) < 0
    stab.append(stable)
    determ.append(det)

mangsV = []
peatsV = []
saltsV = []
for j in range(n):
    mangsV.append(np.real(eigsV[j][0]))
    peatsV.append(np.real(eigsV[j][1]))
    saltsV.append(np.real(eigsV[j][2]))



p1 = plt.figure()
#eigs1 = [np.real(a) for (a,b,c) in eigs]
#eigs2 = [np.real(b) for (a,b,c) in eigs]
#eigs3 = [np.real(c) for (a,b,c) in eigs]

#plt.plot(eigs1)
plt.plot(eigs)

plt.title('Maximum eigenvalue from default to drought & low sea-level')
plt.ylabel('Real part of eigenvalue')
plt.show()

p2 = plt.figure()
plt.plot(mangsV,label='mangroves')
plt.plot(peatsV,label='peats')
plt.plot(saltsV,label='salinity')
plt.legend(loc='best')
plt.title('Eigenvector components from default to drought & low-sea level')
        

#p1 = plt.scatter(x,y)
#plt.xlabel('HydP')
#plt.ylabel('Re(max eigenvalue)')
#plt.show()

# Compute correlations


