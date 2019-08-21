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
#----------------------------------------
# Lotka Volterra parameters
a = 2. # Prey growth rate
b = 0.2 #Background prey mortality rate
c = 1.3 #Bockground predator mortality rate
d = 0.65 #Predation conversion efficiency

# X[0] = prey, x[1] = predator
#----------------------------------------

# Lotka-Volterra equations and Jacobian
#----------------------------------------
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
#----------------------------------------

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

#----------------------------------------
# Plotting phase plane and trajectories
#----------------------------------------

values = np.linspace(0.3, 0.9, 5) # X0 in values between stable points
vcolours = p.cm.autumn_r(np.linspace(0.3, 1, len(values))) #colours for trajectories

f2 = p.figure()
#Plot trajectories
for v, col in zip(values, vcolours):
    X0 = v * X_f1 #start point
    X = integrate.odeint(dX_dt, X0, t) #integrate
    p.plot(X[:,0], X[:,1], lw=3.5*v, color=col, label='X0=(%.f, %.f)' % (X0[0], X0[1]))

#Define grid and compute direction at each ponit
    # Get axes limits
ymax = p.ylim(ymin=0)[1] 
xmax = p.xlim(xmin=0)[1]
nb_points = 20

x = np.linspace(0, xmax, nb_points)
y = np.linspace(0, ymax, nb_points)

X1, Y1 = np.meshgrid(x, y) # create grid
DX1, DY1 = dX_dt([X1,Y1]) # compute infinitesimal growth rate at each point in grid
M = (np.hypot(DX1,DY1) ) #slope for each point
M[M==0] = 1. # Avoid zeros
DX1 /= M #Normalize arrow components
DY1 /= M


#----------------------------------------
# Draw direction fields, use plt's quiver function
# Colours correspond to speed, can also remove normalization

p.title('Phase space of Lotka-Volterra equations')
Q = p.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=p.cm.jet)
p.xlabel('Population of Prey')
p.ylabel('Population of Predators')
p.legend()
p.grid()
p.xlim(0, xmax)
p.ylim(0, ymax)


