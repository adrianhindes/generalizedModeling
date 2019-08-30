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
    
    betaP = data['betaP']
    betaD = data['betaD']
    betaL = data['betaL']
    
    betaA = data['betaA']
    betaR = data['betaR']
    betaE = data['betaE']
    
    hydP = data['hydP']
    propM = data['propM']
    propS = data['propS']
    growM = data['growM']
    
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
    inS = data['inS']
    inM = data['inM']
    outS = data['outS']
    
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

    jac0 =  np.array([R1,R2,R3])
    jac = np.matmul(alphas, jac0)
    
    return jac