# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 23:35:08 2019

@author: Excalibur
Default values and ranges
for generalised modelling parameters
"""
#Mangrove Gain
betaG0 = 0.5
betaP0 = 1-betaG0

# Mangrove loss
betaD0 = 0.3
betaS0 = 0.5
betaL0 = 1-betaD0 - betaS0

# Peat gain
betaA0 = 0.6
betaR0 = 0.1
betaV0 = 1 - betaA0 - betaR0

#Peat loss
betaE0 = 0.5
betaSB0 = 1 - betaE0

# ----------------------
# Elasticity parameters
# ----------------------

# Mangroves
propM0 = 1 
propS0 = -1
growM0 = 1

drownHyd0 = 2
drownM0 = 1 

stressM0 = 1
stressS0 = 2

littM0 = 1.5

# Peat soils
accSed0 = 1 
sedHyd0 = 2
accM0 = 1

retLitt0 = 1
retHyd0 = 1

volGrow0 = 1
volP0 = 1

eroM0 = 1

subsM0 = 1
subsHyd0 = 1
subsP0 = 0.5

# Salinity
inM0 = 1.0 
inS0 = 0.5

outS0 = 1

hydP0 = -1