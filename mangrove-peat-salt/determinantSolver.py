# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 23:43:02 2019

@author: Excalibur
"""

import numpy as np
from numpy import random
from numpy import linalg as LA
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


from parameterDefaults import defaults

sp.init_printing() # Make symbolic expressions look nice

# Bifurcation surface parameters
# Recall: betaG, betaL, betaV and betaSB are defined as 1 - other Betas


points = 1000

# Set ranges


# Helper function for excluding variables from substituting

def subit(expr,defaults,excl):

    repl2 = [(a,b) for a, b in defaults if not any((a==c) for c in excl)]
    
    new = expr.subs(repl2)
    return new

params = []
# Timescale parameters
alphaM, alphaP, alphaS = sp.symbols('alpha_m alpha_p alpha_s')

# Mangrove Betas

betasMang = sp.symbols('beta_p, beta_g, beta_d, beta_s, beta_l')
betaP, betaG, betaD, betaS, betaL = betasMang
params += list(betasMang)
#Peat Betas

betasPeat = sp.symbols('beta_a, beta_r, beta_v, beta_e, beta_sb')
betaA, betaR, betaV, betaE, betaSB = betasPeat
params += list(betasPeat)

# Mangrove Elasticity
elasMang = sp.symbols('grow_m, prop_m, prop_s, drown_hyd, drown_m, stress_m, stress_s, litt_m')
growM, propM, propS, drownHyd, drownM, stressM, stressS, littM = elasMang
params += list(elasMang)

# Peat soils
elasPeat = sp.symbols('acc_sed, sed_hyd, acc_m, ret_litt, ret_hyd, vol_grow, vol_p, ero_m, subs_m, subs_hyd, subs_p')
accSed, sedHyd, accM, retLitt, retHyd, volGrow, volP, eroM, subsM, subsHyd, subsP = elasPeat
params += list(elasPeat)

# Salinity
elasSalt = sp.symbols('in_m, in_s, out_s, hyd_p')
inM, inS, outS, hydP = elasSalt
params += elasSalt

def chSymtoLabel(sym):
    x = str(sym)
    bits = x.split('_')
    if bits[1] == 'sb':
        label = bits[0]+'SB'
    else:
        label = bits[0]+bits[1].title()
    return label

    

p1 = [betaP0, betaG0, betaD0, betaS0, betaL0]
p2 = [betaA0, betaR0, betaV0, betaE0, betaSB0]
p3 = [growM0, propM0, propS0, drownHyd0, drownM0, stressM0, stressS0, littM0]
p4 = [accSed0, sedHyd0, accM0, retLitt0, retHyd0, volGrow0, volP0, eroM0, subsM0, subsHyd0, subsP0]
p5 = [inM0, inS0, outS0, hydP0]

ps = p1+p2+p3+p4+p5

# Tuple list (symbol,default value)
symDefaults = [(sym,defaults[chSymtoLabel(sym)]) for sym in params]

#########
#Bifurcation surface parameters
X = betaA
Y = inS
Z = outS


xMin = 0
xMax = 1

yMin = 0.5
yMax = 1

zAxMin = 0.5
zAxMax = 1


#########
# Jacobian components

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

# Define matrices

alphas = sp.Matrix([[alphaM, 0, 0], [0, alphaP, 0], [0, 0, alphaS]])
jac = sp.Matrix([[dmdm, dmdp, dmds], [dpdm, dpdp, dpds], [dsdm, dsdp, dsds]])
det = jac.det()

# Beta simplifications
# Require additional assumptions if choosing one of trio betas
betaL0 = 0.2
betaS0 = 0.4

betaV0 = 0.2
betaR0 = 0.2

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
if X in betaSet:
    det = fixBetas(det,X)
if Y in betaSet:
    det = fixBetas(det,Y)
    

#det.subs(betaG, subG)
#det.subs(betaL,subL)
#det.subs(betaV,subV)
#det.subs(betaSB,subSB)

saddle = sp.Eq(det,0)
saddleManifold = subit(saddle, symDefaults, [X,Y,Z])

# Surface equation for given X,Y,Z
saddleFunc = sp.solve(saddleManifold, Z)[0]

saddleFun = sp.lambdify((X,Y), saddleFunc)
x = np.linspace(xMin,xMax,points)
y = np.linspace(yMin,yMax,points)

xx, yy = np.meshgrid(x,y)

zz = saddleFun(xx,yy)



fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

ax2.plot_surface(xx, yy, zz)
ax2.set_xlabel('r'+chSymtoLabel(X))
ax2.set_xlim(xMin,xMax)
ax2.set_ylabel('r'+chSymtoLabel(Y))
ax2.set_ylim(yMin,yMax)
ax2.set_zlabel('r'+chSymtoLabel(Z))
#ax2.set_zlim(zAxMin,zAxMax)
plt.show()

from jacobian import computeJac

def checkStability(parX,parY,parZ):
    xSym, x = parX
    ySym, y = parY
    zSym, z = parZ
    # Check stability of system at a point
    labX = chSymtoLabel(xSym)
    labY = chSymtoLabel(ySym)
    labZ = chSymtoLabel(zSym)
    
    data = defaults
    data[labX] = x
    data[labY] = y
    data[labZ] = z
    
    jac = computeJac(data)
    w,v = LA.eig(jac)
    if np.real(max(w)) < 0:
        return 'stable'
    else: return 'unstable'

#dx = 0.05*(xMax-xMin)
#dy = 0.05*(yMax-yMin)

#xt = 0.5*(xMax-xMin)
#yt = 0.5*(yMax-yMin)

#zt = saddleFun(xt,yt)

#res = checkStability((X,xt+dx),(Y,yt+dy),(Z,zt+dx))
#res2 = checkStability((X,xt+dx),(Y,yt+dy),(Z,zt-dx))
#xs = [xt+dx,xt+dx]
#ys = [yt+dy,yt+dy]
#zs = [zt+dx,zt-dx]

#print((res+' at point'+'('+str(xt+dx)+','+str(yt+dy)+','+str(zt+dx)+')'))
#print((res2+' at point'+'('+str(xt+dx)+','+str(yt+dy)+','+str(zt-dx)+')')) 

#ax2.scatter(xs,ys,zs)


#plt.show()
   



