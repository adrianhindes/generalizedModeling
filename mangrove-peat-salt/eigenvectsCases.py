# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from tqdm import tqdm

from jacobianSalt import computeJac
from parameterDefaults import defaults as defaults
from droughtParams import defaults as droughts


# Linspace range
n = 10000

paramLabels = {'alphaM','alphaP','alphaS',
        'betaG','betaP','betaD','betaL','betaS',
        'betaA','betaR','betaV','betaE','betaSB',
        'hydP','propM','propS',\
        'growM','growS',\
        'propPrecip','growPrecip',\
        'drownHyd','drownM','stressM',\
        'stressS','littM','accSed',\
        'sedHyd','accM','retLitt','retHyd',\
        'volGrow','volP','volHyd','volPrecip','eroM',\
        'subsMort','subsHyd','subsP','concS',\
        'concEvapt','concHyd','evaptM','evaptS',\
        'decrS','decrPrecip','precipBeta'}

dataEndpoints = {key:(defaults[key],droughts[key]) for key in paramLabels}
data = {k:np.linspace(v[0],v[1],n) for (k,v) in dataEndpoints.items()}



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


eigs = [] #eigenvalue triplets
eigsV = [] #eigenvector triplets
stab = [] #stability (0,1)
determ = [] # determinant of Jacobian

mangsV = []
peatsV = []
saltsV = []

for j in tqdm(range(n)):
    # Set j-th data point dictionary
    dataJ = {k:v[j] for (k,v) in data.items()}
    jac = computeJac(dataJ)
    
    w, v = LA.eig(jac)
    det = LA.det(jac)
    
    eigMax = np.max(w)
    eigVMax = v[: , w.argmax()]
    
    eigs.append(eigMax)
    eigsV.append(eigVMax)
    stable = np.real(eigMax) < 0
    stab.append(stable)
    determ.append(det)
    
    mangsV.append(np.real(eigsV[j][0]))
    peatsV.append(np.real(eigsV[j][1]))
    saltsV.append(np.real(eigsV[j][2]))


    


# Plot change in maximum eigenvalue



from modelFullvary2 import maxStable, minStable

scenarioParams = set(['betaG','betaD','betaS','betaA','betaR',
                  'betaE','hydP','drownHyd','drownM'])
    
maxStable = {k:v for (k,v) in maxStable.items() if k not in scenarioParams}
minStable = {k:v for (k,v) in minStable.items() if k not in scenarioParams}

eigsStable = []
mangsVstable = []
peatsVstable = []
saltsVstable = []

eigsUnstable = []
mangsVUnstable = []
peatsVUnstable = []
saltsVUnstable = []

for j in tqdm(range(n)):
    # Set j-th data point dictionary
    dataJ = {k:v[j] for (k,v) in data.items()}
    
    for (param,val) in maxStable.items():
        dataJ[param] = val
    
    jac = computeJac(dataJ)
    
    w, v = LA.eig(jac)
    det = LA.det(jac)
    
    eigMax = np.max(w)
    eigVMax = v[: , w.argmax()]
    
    eigsStable.append(np.real(eigMax))
    mangsVstable.append(np.real(eigVMax[0]))
    peatsVstable.append(np.real(eigVMax[1]))
    saltsVstable.append(np.real(eigVMax[2]))
    
    for (param,val) in minStable.items():
        dataJ[param] = val
    
    jac = computeJac(dataJ)
    
    w, v = LA.eig(jac)
    det = LA.det(jac)
    
    eigMax = np.max(w)
    eigVMax = v[: , w.argmax()]
    
    eigsUnstable.append(np.real(eigMax))
    mangsVUnstable.append(np.real(eigVMax[0]))
    peatsVUnstable.append(np.real(eigVMax[1]))
    saltsVUnstable.append(np.real(eigVMax[2]))

p1 = plt.figure()
plt.plot(eigs, label = 'Scenario')
plt.plot(eigsStable, label = 'Scenario w stable parameters')
plt.plot(eigsUnstable, label = 'Scenario w unstable parameters')
plt.legend(loc='best')
#plt.legend(loc='best')

plt.title('Maximum eigenvalue from default to drought & low sea-level')
plt.ylabel(r'Real part of $\lambda$')
plt.xlabel('Stable to ENSO conditions')
plt.show()

# Plot eigenvectors
p2 = plt.figure()
plt.plot(mangsV, label='Mangroves')
#plt.plot(mangsVstable, label='Mangroves Stable')
#plt.plot(mangsVUnstable, label='Mangroves Unstable')
#plt.legend(loc='best')
plt.title('Mangrove Eigenvector component')
plt.show()

p3 = plt.figure()
plt.plot(peatsV,label='Peats')
#plt.plot(peatsVstable, label='Peats Stable')
#plt.plot(peatsVUnstable, label='Peats Unstable')
plt.legend(loc='best')
plt.title('Peats Eigenvector component')
plt.show()

p4 = plt.figure()
plt.plot(saltsV,label='Salinity')
#plt.plot(saltsVstable, label='Salinity Stable')
#plt.plot(saltsVUnstable, label='Salinity Unstable')
plt.legend(loc='best')
plt.title('Salinity Eigenvector component')
plt.show()

