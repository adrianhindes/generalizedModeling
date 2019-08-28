# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 16:47:24 2019

@author: hindesa
"""
from sympy import *
import numpy as np
init_printing()

x, y = symbols('x y')
r, a, b, g, m, k = symbols('r, a, b, gamma, m, k')

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