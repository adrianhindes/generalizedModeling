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
propS = -41
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
data = {'betaP':betaP,'hydP':hydP,
        'drownHyd':drownHyd,'stressS':stressS,
        'littM':littM,'sedHyd':sedHyd,'retHyd':retHyd,
        'subsHyd':subsHyd,'subsP':subsP,'inS':inS,
        'eigs':[], 'eigsV':[],'stab':[]}

for j in tqdm(range(n)):
    
    # Remap parameters
    betaP = data['betaP'][j]
    hydP = data['hydP'][j]
    drownHyd = data['drownHyd'][j]
    stressS = data['stressS'][j]
    littM = data['littM'][j]
    sedHyd = data['sedHyd'][j]
    retHyd = data['retHyd'][j]
    subsHyd = data['subsHyd'][j]
    subsP = data['subsP'][j]
    inS = data['inS'][j]
    
    # Define Jacobian matrix elements
    # Note syntax here is not strictly correct - dmdm == (dm/dt)/dm
    dmdm = betaP*propM +(1-betaP)*growM -betaS*stressM -betaD*drownM -betaL*littM
    dmdp = -betaD*hydP*drownHyd
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
    
    data['eigs'].append(w)
    data['eigsV'].append(v)
    data['stab'].append(stability(w))
    


x = data['hydP']
ys = data['eigs']
y0 = []
y1 = []
y2 = []

yz = []
for j in range(n):
    y0.append(np.real(ys[j][0]))
    y1.append(np.real(ys[j][1]))
    y2.append(np.real(ys[j][2]))
    
    yz.append(np.max(np.real(ys[j])))

p1 = plt.scatter(range(n),yz)

# Compute correlations

        
betaPCor = np.corrcoef(data['betaP'],data['stab'])[0,1]
hydPCor = np.corrcoef(data['hydP'],data['stab'])[0,1]
#propSCor = np.corrcoef(data['propS'],data['stab'])[0,1]
drownHydCor = np.corrcoef(data['drownHyd'],data['stab'])[0,1]
stressSCor = np.corrcoef(data['stressS'],data['stab'])[0,1]
littMCor = np.corrcoef(data['littM'],data['stab'])[0,1]
sedHydCor = np.corrcoef(data['sedHyd'],data['stab'])[0,1]
retHydCor = np.corrcoef(data['retHyd'],data['stab'])[0,1]
subsHydCor = np.corrcoef(data['subsHyd'],data['stab'])[0,1]
subsPCor = np.corrcoef(data['subsP'],data['stab'])[0,1]
inSCor = np.corrcoef(data['inS'],data['stab'])[0,1]

corrs = [betaPCor, hydPCor, drownHydCor, stressSCor, littMCor,
         sedHydCor, retHydCor, subsHydCor, subsPCor, inSCor]
labels = ['betaP','hydP','drownHyd','stressS','littM','sedHyd',
          'retHyd','subsHyd','subsP','inS']

p2 = plt.bar(np.arange(11), corrs)
plt.ylabel('Correlation Coefficient')
plt.xlabel('Parameter')
plt.xticks(np.arange(11), labels,rotation=70)
plt.show()
