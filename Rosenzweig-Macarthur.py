# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 16:47:24 2019

@author: hindesa
"""
from sympy import *
from sympy.plotting import plot3d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
init_printing()

x, y = symbols('x y')
r, a, b, g, m, k = symbols('r, a, b, gamma, m, k')

X = k
Y = a
Z = m

rangeX = (X,-1,1)
rangeY = (Y,-1,1)

v = Matrix([x, y])

xdot = r*x*(1-x/k) - (a*x*y)/(b+x)
ydot = g*(a*x*y)/(b+x) - m*y

xdot0 = Eq(xdot, 0)

ydot0 = Eq(ydot, 0)

steadies = solve((xdot0, ydot0), x, y)

rmac = Matrix([xdot, ydot])
jac = rmac.jacobian(v)

fix1 = steadies[0]
fix2 = steadies[1]
fix3 = steadies[2]

jac1 = jac.subs({x: fix1[0], y: fix1[1]})
jac2 = jac.subs({x: fix2[0], y: fix2[1]})
jac3 = jac.subs({x: fix3[0], y: fix3[1]})

det1 = jac1.det()
det2 = jac2.det()
det3 = jac3.det()

r0 = 1
a0 = 0.5
b0 = 3.0
g0 = 2.2
m0 = 1
k0 = 3
repl = [(r,r0), (a,a0), (b,b0), (g,g0), (m,m0), (k,k0)]

def subit(expr,excl):

    repl2 = [(a,b) for a, b in repl if not any((a==c) for c in excl)]
    
    new = expr.subs(repl2)
    return new

check = {X,Y,Z}
XY = [X,Y]

if check.issubset(det2.atoms()):
    
    surf2 = solve(Eq(det2,0),Z)[0]
    
    surf = subit(surf2, [X,Y])
    
    plot3d(surf,rangeX,rangeY,xlabel=X,ylabel=Y,title=Z)
    
