# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 23:35:08 2019

@author: Excalibur
Default (healthy) values and ranges
for generalised modelling parameters
"""
alphaM0 = 1/6.
alphaP0 = 1/3.
alphaS0 = 6.

alphas = {'alphaM':alphaM0, 'alphaP':alphaP0, 'alphaS':alphaS0}
#Mangrove Gain
betaG0 = 0.5
betaP0 = 1-betaG0

# Mangrove loss
betaD0 = 0.4
betaS0 = 0.4
betaL0 = 1-betaD0 - betaS0

# Peat gain
betaA0 = 0.6
betaR0 = 0.1
betaV0 = 1 - betaA0 - betaR0

#Peat loss
betaE0 = 0.8
betaSB0 = 1 - betaE0

betas = {'betaG':betaG0, 'betaP':betaP0, 'betaD':betaD0, 'betaS':betaS0, 'betaL':betaL0,
         'betaA':betaA0, 'betaR':betaR0, 'betaV':betaV0, 'betaE':betaE0, 'betaSB':betaSB0}

# ----------------------
# Elasticity parameters
# ----------------------

# Mangroves
propM0 = 1 
propS0 = -2
growM0 = 1
growS0 = -1

propPrecip = 2
growPrecip = 1
evaptM = 0.5
precipBeta = 0.5

drownHyd0 = 2
drownM0 = 1 

stressM0 = 1
stressS0 = 2

littM0 = 1.5

mangs = {'propM':propM0, 'propS':propS0, 'growM':growM0,'growS':growS0, 'drownHyd':drownHyd0, \
         'drownM':drownM0,'stressM':stressM0, 'stressS':stressS0, 'littM':littM0,\
         'propPrecip':propPrecip,'growPrecip':growPrecip,\
         'evaptM':evaptM,'precipBeta':precipBeta}



# Peat soils
accSed0 = 1
sedHyd0 = 2
accM0 = 1.5

retLitt0 = 1
retHyd0 = -1

volGrow0 = 1
volP0 = 1
volHyd0 = 1
volPrecip = 0.5

eroM0 = -1

subsMort0 = 1.5
subsHyd0 = 1
subsP0 = 0.5

hydP0 = -1

peats = {'accSed':accSed0, 'sedHyd':sedHyd0, 'accM':accM0,\
         'retLitt':retLitt0, 'retHyd':retHyd0, 'volGrow':volGrow0,
         'volP':volP0,'volPrecip':volPrecip, 'eroM':eroM0, 'subsMort':subsMort0,\
         'subsHyd':subsHyd0, 'subsP':subsP0, 'hydP':hydP0,'volHyd':volHyd0}

# Salinity
concEvapt = 1.0 
concS = 0.8
concHyd = 1

decrS = 0.5
decrPrecip = 1

evaptS = -1
evaptM = 1.5

salts = {'concEvapt':concEvapt,'concHyd':concHyd, 'concS':concS, 'decrS':decrS,
         'decrPrecip':decrPrecip,'evaptS':evaptS}

defaults = {**alphas, **betas, **mangs, **peats, **salts}