# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 12:52:44 2019

@author: Excalibur
"""

import numpy as np
import numpy.linalg as LA
from numpy import random
import matplotlib.pyplot as plt
from parameterDefaults import defaults
from jacobianSalt import computeJac
from parameterRanges import ranges
from tqdm import tqdm


dp = 0.01
param = 'betaA'

J1 = computeJac(defaults)
w1, v = LA.eig(J1)

#perturb = defaults
#perturb[param] = defaults[param] + dp

#J2 = computeJac(perturb) 
#w2, v = LA.eig(J2)

#dw2 = w2-w1
#fracPar = defaults[param]/dp

#eigSens = [fracPar*dwi for dwi in dw2]

mangs = {'propM', 'propS', 'growM','growS', 'drownHyd','drownM',
         'stressM', 'stressS', 'littM','propPrecip','growPrecip',
         'evaptM','precipEvapt'}

peats = {'accSed', 'sedHyd', 'accM','retLitt', 'retHyd', 'volGrow',
         'volP','volPrecip', 'eroM', 'subsMort', 'subsHyd', 'subsP',
         'hydP','volHyd'}

salts = {'concEvapt','concHyd', 'concS', 'decrS','decrPrecip','evaptS'}

elasPars = mangs.union(peats).union(salts)
nRuns = 1000

elasEigs = {par:np.zeros(nRuns) for par in elasPars}

for par in elasPars:
    r0 = ranges[par][0]
    r1 = ranges[par][1]
    parRange = np.linspace(r0,r1,nRuns)
    
    parSet = defaults
    for j in range(nRuns):
        parVal = parRange[j]
        for p in elasPars:
            parSet[p] = random.uniform(ranges[p][0], ranges[p][1])

        parSet[par] = parVal
        
        J = computeJac(parSet)
        w, v = LA.eig(J)
        
        maxW = np.max(np.real(w))
        elasEigs[par][j] = maxW
        
p1 = plt.figure()
for par in mangs:
    plt.plot(range(nRuns), elasEigs[par], label=par, marker='+')
    
plt.legend(loc='best')

p2 = plt.figure()
for par in peats:
    plt.plot(range(nRuns), elasEigs[par], label=par, marker='+')
    
plt.legend(loc='best')

p3 = plt.figure()
for par in salts:
    plt.plot(range(nRuns), elasEigs[par], label=par, marker='+')
    
plt.legend(loc='best')     