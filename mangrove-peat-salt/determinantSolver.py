# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 23:43:02 2019

@author: Excalibur
"""

import numpy as np
from numpy import random
from numpy import linalg as LA
import sympy as sp
from sympy.printing import latex
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

from parameterDefaults import defaults
from droughtParams import defaults as drought
from parameterRanges import ranges
from jacobianSalt import computeJac
#sp.init_printing() # Make symbolic expressions look nice

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
alphas = sp.symbols('alpha_m alpha_p alpha_s')
alphaM, alphaP, alphaS = alphas
params += list(alphas)

# Mangrove Betas

betasMang = sp.symbols('beta_p, beta_g, beta_d, beta_s, beta_l')
betaP, betaG, betaD, betaS, betaL = betasMang
params += list(betasMang)
#Peat Betas

betasPeat = sp.symbols('beta_a, beta_r, beta_v, beta_e, beta_sb')
betaA, betaR, betaV, betaE, betaSB = betasPeat
params += list(betasPeat)

# Mangrove Elasticity
elasMang = sp.symbols('grow_m, grow_s, prop_m, prop_s, drown_hyd, drown_m,\
                      stress_m, stress_s,litt_m, \
                      prop_precip, grow_precip, precip_beta')
growM, growS, propM, propS, drownHyd, drownM, stressM, stressS, littM,\
                     propPrecip, growPrecip, precipBeta = elasMang
params += list(elasMang)

# Peat soils
elasPeat = sp.symbols('acc_sed, sed_hyd, acc_m, ret_litt, ret_hyd, vol_grow,\
                      vol_p, ero_m, subs_mort, subs_hyd, subs_p, vol_hyd, vol_precip')
accSed, sedHyd, accM, retLitt, retHyd, volGrow, volP, eroM, subsMort, subsHyd, subsP,volHyd, volPrecip = elasPeat
params += list(elasPeat)

# Salinity
elasSalt = sp.symbols('conc_evapt, conc_hyd, evapt_m, decr_precip, conc_s, decr_s, hyd_p, evapt_s')
concEvapt, concHyd, evaptM, decrPrecip,concS, decrS, hydP, evaptS = elasSalt
params += elasSalt

def chSymtoLabel(sym):
    x = str(sym)
    bits = x.split('_')
    if bits[1] == 'sb':
        label = bits[0]+'SB'
    else:
        label = bits[0]+bits[1].title()
    return label

    

#p1 = [betaP0, betaG0, betaD0, betaS0, betaL0]
#p2 = [betaA0, betaR0, betaV0, betaE0, betaSB0]
#p3 = [growM0, propM0, propS0, drownHyd0, drownM0, stressM0, stressS0, littM0]
#p4 = [accSed0, sedHyd0, accM0, retLitt0, retHyd0, volGrow0, volP0, eroM0, subsM0, subsHyd0, subsP0]
#p5 = [inM0, inS0, outS0, hydP0]

#ps = p1+p2+p3+p4+p5

# Tuple list (symbol,default value)
symDefaults = [(sym,defaults[chSymtoLabel(sym)]) for sym in params]
symDroughts = [(sym,drought[chSymtoLabel(sym)]) for sym in params]


############
#----------#
############

#Bifurcation surface parameters
# Top correlation parameters to check
# hydP, decrS, betaA, betaSB, concEvapt, evaptM
# concS, betaE, betaV


X = betaS
Y = stressS
Z = concHyd

blacklist = [betaP, betaL, betaR] 
check = [(par in blacklist) for par in [X,Y,Z]]
if any(check): sys.exit("One or more parameters is an already defined beta")

showDrought = False
############
#----------#
############

xMin = ranges[chSymtoLabel(X)][0]
xMax = ranges[chSymtoLabel(X)][1]

yMin = ranges[chSymtoLabel(Y)][0]
yMax = ranges[chSymtoLabel(Y)][1]-0.1

zAxMin = ranges[chSymtoLabel(Z)][0]
zAxMax = ranges[chSymtoLabel(Z)][1]

truncate = True
#########
# Jacobian components
mortD = betaD/(betaD+betaS)
mortS = betaS/(betaD+betaS)

dPropdM = propM+propPrecip*precipBeta*evaptM
dGrowdM = growM+growPrecip*precipBeta*evaptM

# Beta substitutions
betaP = 1 - betaG
betaL = 1 - betaD - betaS

betaR = 1 - betaA - betaV
betaE = 1 - betaSB


dmdm = betaP*dPropdM +betaG*dGrowdM-betaS*stressM -betaD*drownM -betaL*littM
         
dmdp = -1*betaD*hydP*drownHyd

dPropdS = propS+propPrecip*precipBeta*evaptS
dGrowdS = growS+growPrecip*precipBeta*evaptS

dmds = betaP*dPropdS + betaG*dGrowdS-betaS*stressS

dVoldM = volGrow*(growM+growPrecip*precipBeta*evaptM)+volPrecip*precipBeta*evaptM
dSubsdM = subsMort*(mortD*drownM+mortS*stressM)

dpdm = betaA*accM +betaR*retLitt*littM + betaV*dVoldM - betaE*eroM -betaSB*dSubsdM

dVoldP = volHyd*hydP+volP
dSubsdP = subsHyd*hydP + subsP
            
dpdp = hydP*(betaA*accSed*sedHyd + betaR*retHyd)+betaV*dVoldP-betaSB*dSubsdP


dVoldS = volGrow*(growS+growPrecip*precipBeta*evaptS)+volPrecip*precipBeta*evaptS
dSubsdS = subsMort*(mortS*stressS)

dpds = betaV*dVoldS - betaSB*dSubsdS
    
dsdm = evaptM*(concEvapt - decrPrecip*precipBeta)
dsdp = concHyd*hydP
dsds = concEvapt*evaptS - decrPrecip*precipBeta*evaptS

# Define matrices

alphas = sp.Matrix([[alphaM, 0, 0], [0, alphaP, 0], [0, 0, alphaS]])
jac = sp.Matrix([[dmdm, dmdp, dmds], [dpdm, dpdp, dpds], [dsdm, dsdp, dsds]])
jac2 = alphas*jac
det = jac2.det()



saddle = sp.Eq(det,0)
saddleManifold = subit(saddle, symDefaults, [X,Y,Z])





# Surface equation for given X,Y,Z
saddleFunc = sp.solve(saddleManifold, Z)[0]

print(str(Z) +'=' + str(saddleFunc))

saddleFun = sp.lambdify((X,Y), saddleFunc)

xs = np.linspace(xMin,xMax,points)
ys = np.linspace(yMin,yMax,points)

xx, yy = np.meshgrid(xs,ys)
zz = saddleFun(xx,yy)

if showDrought: 
    saddleManifold2 = subit(saddle, symDroughts, [X,Y,Z])
    saddleFunc2 = sp.solve(saddleManifold2, Z)[0]
    print(str(Z) +'=' + str(saddleFunc2))
    saddleFun2 = sp.lambdify((X,Y), saddleFunc2)
    zz2 = saddleFun2(xx,yy)



fig2=plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

if truncate == True:
    for i in range(len(xx)):
        for j in range(len(yy)):
            if (zz[j,i] < zAxMin) or (zz[j,i] > zAxMax): zz[j,i] = np.nan
            if showDrought:
                if (zz2[j,i] < zAxMin) or (zz2[j,i] > zAxMax): zz2[j,i] = np.nan
            

        
s1 = ax2.plot_surface(xx, yy, zz, alpha=0.9, edgecolor='none')

s1._facecolors2d=s1._facecolors3d
s1._edgecolors2d=s1._edgecolors3d

if showDrought == True:
    s2 = ax2.plot_surface(xx, yy, zz2, label='Drought', alpha=0.7, edgecolor='none')
    s2._facecolors2d=s2._facecolors3d
    s2._edgecolors2d=s2._edgecolors3d
    ax2.legend(loc = 'upper right')
    
ax2.set_xlabel(r'$'+latex(X)+'$')
ax2.set_xlim(xMin,xMax)
ax2.set_ylabel(r'$'+latex(Y)+'$')
ax2.set_ylim(yMin,yMax)
ax2.set_zlabel(r'$'+latex(Z)+'$')
#plt.title(r'Bifurcation Surface of $('+latex(X)+','+latex(Y)+','+latex(Z)+')$')

# Plot stable parameter triplet from defaults (assuming defaults are stable)
stX = defaults[chSymtoLabel(X)]
stY = defaults[chSymtoLabel(Y)]
stZ = defaults[chSymtoLabel(Z)]

ax2.scatter(stX,stY,stZ, color='red', label = 'Stable config')
#ax2.legend(loc='lower left')

if truncate == True:
    ax2.set_zlim(zAxMin,zAxMax)




plt.show()




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


gradNorm = [-saddleFunc.diff(X),-saddleFunc.diff(Y),1]
gradNorm = np.multiply(gradNorm, 1/10)
x1 = random.uniform(xMin,xMax)
y1 = random.uniform(yMin,yMax)

gradNorm[0] = gradNorm[0].subs([(X,x1),(Y,y1)])
gradNorm[1] = gradNorm[1].subs([(X,x1),(Y,y1)])

revGradNorm = np.multiply(gradNorm, -1)


p1 = (x1,y1,saddleFun(x1,y1))
p2 = [a+b for a,b in zip(gradNorm,p1)]
p3 = [a+b for a,b in zip(revGradNorm,p1)]

#res1 = checkStability((X,p2[0]),(Y,p2[1]),(Z,p2[2]))
#print(res1 + " above" )
#res2 = checkStability((X,p3[0]),(Y,p2[1]),(Z,p3[2]))
#print(res2 + " below")

#res = checkStability((X,xt+dx),(Y,yt+dy),(Z,zt+dx))
#res2 = checkStability((X,xt+dx),(Y,yt+dy),(Z,zt-dx))

#print((res+' at point'+'('+str(xt+dx)+','+str(yt+dy)+','+str(zt+dx)+')'))
#print((res2+' at point'+'('+str(xt+dx)+','+str(yt+dy)+','+str(zt-dx)+')')) 

#ax2.scatter(xs,ys,zs,color='red')


#plt.show()
   



