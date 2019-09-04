# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 13:17:03 2019

@author: hindesa
return Jacobian for mangroves model given parameter set
"""
import numpy as np

def computeJac(data):
    alphaM = data['alphaM']
    alphaP = data['alphaP']
    alphaS = data['alphaS']
    
    betaG = data['betaG']
    betaP = data['betaP']
    
    betaP = data['betaP']
    betaD = data['betaD']
    betaL = data['betaL']
    
    betaA = data['betaA']
    betaR = data['betaR']
    betaE = data['betaE']
    
    hydP = data['hydP']
    propM = data['propM']
    propS = data['propS']
    growS = data['growS']
    growM = data['growM']
    
    propPrecip = data['propPrecip']
    growPrecip = data['growPrecip']
    
    
    drownHyd = data['drownHyd']
    drownM = data['drownM']
    stressM = data['stressM']
    stressS = data['stressS']
    
    littM = data['littM']
    accSed = data['accSed']
    sedHyd = data['sedHyd']
    accM = data['accM']
    retLitt = data['retLitt']
    retHyd = data['retHyd']
    volGrow = data['volGrow']
    volP = data['volP']
    
    eroM = data['eroM']
    subsM = data['subsM']
    subsHyd = data['subsHyd']
    subsP = data['subsP']
    
    
    concS = data['concS']
    concEvapt = data['concEvapt']
    evaptM = data['evaptM']
    concHyd = data['concHyd']
    
    decrS = data['decrS']
    decrPrecip = data['decrPrecip']
    precipEvapt = data['precipEvapt']
    
    
    
    # Define Jacobian matrix elements
    # Note syntax here is not strictly correct - dmdm == (dm/dt)/dm
    
    betaG = 1-betaP
    betaS = 1-betaD-betaL
    betaSB = 1-betaE
    betaV = 1-betaA-betaR
    
    precipM = evaptM*precipEvapt
    
    dmdm = betaP*(propM + propPrecip*precipM) +betaG*(growM+growPrecip*precipM)\
            +betaS*stressM -betaD*drownM -betaL*littM
    dmdp = -1*betaD*hydP*drownHyd
    dmds = betaP*propS +betaG*growS - betaS*stressS

    dpdm = betaA*accM +betaR*retLitt*littM + betaV*volGrow*growM\
            -betaE*eroM -betaSB*subsM
        
    dpdp = hydP*(betaA*accSed*sedHyd + betaR*retHyd-betaSB*subsHyd)\
            +betaV*volP-betaSB*subsP
    dpds = 0
    
    dsdm = concEvapt*evaptM - decrPrecip*precipEvapt*evaptM
    dsdp = concHyd*hydP
    dsds = concS - decrS
    
    # alpha paramater array
    alphas = np.array([ [alphaM, 0, 0], [0, alphaP, 0], [0, 0, alphaS]])
    alphas = alphas.astype(float)

    R1 = [dmdm, dmdp, dmds]
    R2 = [dpdm, dpdp, dpds]
    R3 = [dsdm, dsdp, dsds]

    jac0 =  np.array([R1,R2,R3])
    jac0 = jac0.astype(float)
    jac = np.matmul(alphas, jac0)
    
    return jac