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

# M = mangrove biomass (kg), P = peat soil elevation (mm)
# S = soil salinity (ppm?)

# Linspace range
n = 10000

# ---------------
# Timescales
# ---------------
# 1/years turnover rates
# total guesses at the moment
alphaM = 1/6
alphaP = 1/3.
alphaS = 3.



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
hydP = random.uniform(-2.0, 2.0, n)
# Mangroves
propM = random.uniform(1, 2, n) 
propS = random.uniform(-1, 0.0, n)
growM = random.uniform(1, 2, n)

drownHyd = random.uniform(0.0, 5.0, n)
drownM = random.uniform(0, 1, n) 

stressM = random.uniform(0, 2, n)
stressS = random.uniform(0.0, 2.0, n)

littM = random.uniform(1.0, 2.0, n)

# Peat soils
accSed = random.uniform(1, 2, n) 
sedHyd = random.uniform(0.0, 2.0, n)
accM = random.uniform(1, 2, n)

retLitt = random.uniform(0, 2, n)
retHyd = random.uniform(-2.0, 0.0, n)

volGrow = random.uniform(1, 2, n)
volP = random.uniform(0.5, 1.5, n)

eroM = random.uniform(1, 2, n)

subsM = random.uniform(1, 3, n)
subsHyd = random.uniform(0.0, 1.0, n)
subsP = random.uniform(0.5, 1.5, n)

# Salinity
inM = random.uniform(0,2,n)
inS = random.uniform(0.0, 1.0, n)

outS = 1

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
        'subsHyd':subsHyd,'subsP':subsP,'inS':inS,'inM':inM}

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

eigVDirsStable = {'+++':0,'++-':0,'+-+':0,'+--':0,'-++':0,'-+-':0,'--+':0,'---':0}
eigVDirsUnstable = {'+++':0,'++-':0,'+-+':0,'+--':0,'-++':0,'-+-':0,'--+':0,'---':0}

def plusMin(val):
    if val < 0:
        ret = '-'
    elif val > 0:
        ret = '+'
    else: ret = '0'
    return ret

stableIndices = []
unstableIndices = []

for j in range(n):
    if stab[j]:
        stableIndices.append(j)
for j in range(n):
    if not(stab[j]):
        unstableIndices.append(j)



for j in range(len(stableIndices)):
    i = stableIndices[j]
    key = ''
    key += plusMin(mangsV[i])
    key += plusMin(peatsV[i])
    key += plusMin(saltsV[i])
    eigVDirsStable[key] += 1

for j in range(len(unstableIndices)):
    i = unstableIndices[j]
    key = ''
    key += plusMin(mangsV[i])
    key += plusMin(peatsV[i])
    key += plusMin(saltsV[i])
    eigVDirsUnstable[key] += 1



        

#p1 = plt.scatter(x,y)
#plt.xlabel('HydP')
#plt.ylabel('Re(max eigenvalue)')
#plt.show()

# Compute correlations


