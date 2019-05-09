# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:47:34 2019

@author: hindesa

Lotka-Volterra specific case
Following SciPy Tutorial:
    https://scipy-cookbook.readthedocs.io/items/LoktaVolterraTutorial.html
"""

import numpy as np
import pylab as p
from scipy import integrate

# Lotka Volterra parameters
a = 1. # Prey growth rate
b = 0.1 #Background prey mortality rate
c = 1.5 #Bockground predator mortality rate
d = 0.75 #Predaction conversion efficiency

# X[0] = prey, x[1] = predator

def dX_dt(A, t=0):
    """ Return the growth rate of prey and predator population """
    x = A[0]
    y = A[1]
    matrix = np.array([a*x - b*x*y,
                    -c*y +d*b*x*y])
    return matrix

# Position equilbrium, in this case obtained algebraically
X_f0 = np.array([0., 0.]) # Extinction
X_f1 = np.array([ c/(d*b), a/b]) #Periodic equilibrium

# Linearize around these fixed points using Jacobian

def jacX(A, t=0):
    """ Return Jacobian matrix evaluated at X """
    x = A[0]
    y = A[1]
    matrix = np.array([[a-b*y, -b*x],
                      [d*b*y, d*b*x - c]])
    return matrix

# Calculate linearizations
A_f0 = jacX(X_f0)
A_f1 = jacX(X_f1)

print("Extinction Jacobian: \n" + str(A_f0))
print("Periodic Jacobian: \n" + str(A_f1))

# Calculate Eigenvalues of periodic solution
lambda1, lambda2 = np.linalg.eigvals(A_f1)
period = 2*np.pi/np.abs(lambda1)
print("Period = " + str(period))

# Integrate system of ODEs
t = np.linspace(0, 15, 1000) #time
x = 10 # Starting no. Prey
y = 5 # Starting no. Predators
X0 = [x,y]

print('Integrating ODE...')
X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
print(infodict['message'])

#Plot Periodic solution
prey, predator = X.T #Assign values by transposing integration result
f1 = p.figure() #Initialize Figure

#Plot data
p.plot(t, prey, 'b-', label='Prey')
p.plot(t, predator, 'r-', label='Predators')
p.grid() #gridlines
p.legend(loc='best') #legend
p.xlabel('Time') #labels
p.ylabel('Population')

p.title('Lotka-Volterra model periodic steady state')

# Importantly, you cannot get a time series out of
# generalised models because there is no functional forms specified
# for an ode solver to do its thing


