# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 13:43:31 2019

@author: Adrian Hindes
Following http://systems-sciences.uni-graz.at/etextbook/sw2/crittrans.html
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import integrate
import sympy as sp
from tqdm import tqdm
n = 1000
fig = plt.figure()

# define model parameters
h = 1.0; b = 1.0; p = 12; a=0.7

def hill(x):
    y = x**p/(x**p+h**p)
    return y


def model(x,r):
    y = a-b*x+r*hill(x)
    return y

def deriv(x,r):
    dgdx = -b+(r*p*x**(p-1)*h**p)/((x**p + h**p)**2)
    return dgdx

rMin = 0
rMax = 3
RS = np.linspace(0,3,n)
XS = np.linspace(0.5,3.5,3000)
tol = 1E-4 # tolerance for close to 0


EX1 = []; EX2 = []; EX3 = []

for r in tqdm(RS):
    for x in XS:
        gx = a-b*x + (r*x**p)/(x**p + h**p)
        # If near equilibrium
        if -tol < gx < tol:
            #calculate g'(x) (= dg(x)/dx) to "check" stability:
            Dgx = -b+(r*p*x**(p-1)*h**p)/((x**p + h**p)**2)
            #write x-coordinate of extremum and depending value of parameter a in an array:
            if Dgx > 0: EX2.append([x,r])      #array for unstable extrema
            if Dgx < 0:                        
                if x < 1: EX1.append([x,r])    #array for stable extrema with x<1
                if x > 1: EX3.append([x,r])    #array for stable extrema with x>1
    #plot f(x) for specific values of a:

EX1 = list(zip(*sorted(EX1,key=lambda x: x[0]))); EX1x = EX1[1]; EX1y = EX1[0]
EX2 = list(zip(*sorted(EX2,key=lambda x: x[0]))); EX2x = EX2[1]; EX2y = EX2[0]
EX3 = list(zip(*sorted(EX3,key=lambda x: x[0]))); EX3x = EX3[1]; EX3y = EX3[0]

def hysteresis(EX1x,EX1y,EX2x,EX2y,EX3x,EX3y,x_label):

    fig = plt.figure(figsize=(12,8))
    mpl.rc('font', size = 14)
    ax = fig.add_subplot(1,1,1)
    
    #plot curve containing the extrema, splitted in lower (solid line), central (dashdotted line)
    #and upper part (solid line) in black color:
    ax.plot(EX1x,EX1y,color = 'k', linewidth=3, label = 'stable equilibria')
    ax.plot(EX2x,EX2y,color = 'k', linewidth=3, linestyle = 'dashdot', label = 'unstable equilibria')
    ax.plot(EX3x,EX3y,color = 'k', linewidth=3)
    
    #plot development of x_extremum starting at maximum value of parameter a
    #find x-coordinate in EX3x at which critical transition occurs:
    for index3, value3 in enumerate(EX3x):
        if -0.005 < (EX1x[-1] - value3) < 0.005:
            break
    line, = ax.plot(EX1x,EX1y,color = 'c', linewidth=8, linestyle='-',
            label='development of $x_{equilibrium}$ for increasing \nvalues of '
            + x_label + ', starting at ' + x_label + '$_{min}$')
    #line.set_dashes((dash_length,dash_distance)) allows to vary the dash-style
    line.set_dashes((1,5))
    line, = ax.plot([EX1x[-1],EX1x[-1]],[EX1y[-1],EX3y[index3]],color = 'c', linewidth=8, linestyle='-')
    line.set_dashes((1,5))
    line, = ax.plot(EX3x[index3:],EX3y[index3:],color = 'c', linewidth=8, linestyle='-')
    line.set_dashes((1,5))
    
    #plot development of x_extremum starting at minimum value of parameter a
    #find x-coordinate in EX1x at which critical transition occurs:
    for index1, value1 in enumerate(EX1x):
        if -0.005 < (EX3x[0] - value1) < 0.005:
            break
    line, = ax.plot(EX3x,EX3y,color = 'r', linewidth=8, linestyle='-',
            label='development of $x_{equilibrium}$ for decreasing \nvalues of '
            + x_label + ', starting at ' + x_label + '$_{max}$')
    line.set_dashes((1,5))
    line, = ax.plot([EX3x[0],EX3x[0]],[EX3y[0],EX1y[index1]],color = 'r', linewidth=8, linestyle='-')
    line.set_dashes((1,5))
    line, = ax.plot(EX1x[:index1],EX1y[:index1],color = 'r', linewidth=8, linestyle='-')
    line.set_dashes((1,5))
    
    return ax, index1, index3

ax, index1, index3 = hysteresis(EX1x,EX1y,EX2x,EX2y,EX3x,EX3y,'$r$')
ax.annotate('', xy=(1.2, 2.1), xytext=(1.5, 2.4),arrowprops=dict(facecolor='r',lw=0),)
ax.annotate('', xy=(1.5, 0.58), xytext=(1.2, 0.565),arrowprops=dict(facecolor='c',lw=0),)
ax.set_title('Hysteresis depicting the development of equilibria (x-cor) depending on variation of r')
ax.grid();
ax.set_xlabel('$r$', fontsize=18);
ax.set_ylabel('$x$', fontsize=18);
ax.set_xlim([EX1x[0],EX3x[-1]])
ax.set_ylim([EX1y[0]-0.3,EX3y[-1]+0.2])
ax.legend(loc='best')

#plt.plot(EX1x,EX1y,color = 'k', linewidth=2, label = 'stable equilibria');
#central part of the curve (unstable extrema):
#plt.plot(EX2x,EX2y,color = 'k', linewidth=2, label = 'unstable equilibria', linestyle = 'dashdot');
#upper part of the curve (stable extrema):
#plt.plot(EX3x,EX3y,color = 'k', linewidth=2);

#configure plot properties:
#plt.title('equilibria of x(t) for specific values of r')
#plt.grid();
#plt.xlabel('$r$', fontsize=18)
#plt.ylabel('$x$', fontsize=18)
#plt.legend(loc='best', bbox_to_anchor=(1, 0.5))