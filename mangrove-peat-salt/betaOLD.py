# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 15:08:28 2019

@author: hindesa
"""

# Beta simplifications
# Require additional assumptions if choosing one of trio betas
betaL0 = 0.8
betaS0 = 0.1

betaV0 = 0.25
betaR0 = 0.05

subG = 1 - betaP
subP = 1 - betaG

subD = 1 - betaL0 - betaS
subS = 1 - betaD - betaL0
subL = 1 - betaD - betaS

subA = 1 - betaR - betaV0
subR = 1 - betaA - betaV0
subV = 1 - betaA - betaR

subE = 1 - betaSB
subSB = 1 - betaE


def fixBetas(expr,beta):
    if beta == betaG:
        return expr.subs(betaP,subP)
    elif beta == betaP:
        return expr.subs(betaG,subG)
    
    elif beta == betaD:
        return expr.subs(betaS,subS)
    elif beta == betaS:
        return expr.subs(betaD,subD)
    elif beta == betaL:
        return expr.subs(betaD,(1-betaL-betaS0))
    
    elif beta == betaA:
        return expr.subs(betaR,subR)
    elif beta == betaR:
        return expr.subs(betaA,subA)
    elif beta == betaV:
        return expr.subs(betaA,(1-betaR0-betaV))
    
    elif beta == betaE:
        return expr.subs(betaSB,subSB)
    elif beta == betaSB:
        return expr.subs(betaE,subE)

betaSet = betasMang + betasPeat

#if X in betaSet: det = fixBetas(det,X)
#if Y in betaSet: det = fixBetas(det,Y)
        
#det.subs(betaG, subG)
#det.subs(betaL,subL)
#det.subs(betaV,subV)
#det.subs(betaSB,subSB)