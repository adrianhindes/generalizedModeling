# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:10:44 2019

@author: Excalibur
Budworm model from Strogatz
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def y1(x,r,k):
    y = r*(1-x/k)
    return y

def y2(x):
    y = x/(1+x**2)
    return y

def potent(x,r,k):
    y = 1/(3*k)*x**3 - (r/2)*x**2 + x - np.arctan(x)
    return y

n = 1000
xs = np.linspace(0,10,n)

c1 = y2(xs)

c2 = y1(xs,0.3,10)
c3 = y1(xs,0.38,10)
c4 = y1(xs,0.5,10)

p = potent(xs,0.4,10)
plt.plot(xs,p)


#plt.plot(xs,c1)
#plt.plot(xs,c2)
#plt.plot(xs,c3)
#plt.plot(xs,c4)

#plt.ylim([0,max(c1)+0.2])
#plt.xlim([0,10])

def test(x,r,k):
    error = y1(x,r,k) - y2(x)
    return error


def r(x):
    y = (2*x**3)/((1+x**2)**2)
    return y

def k(x):
    y = (2*x**3)/(x**2-1)
    return y

xs = np.linspace(1.01,20,n)

rs = r(xs)
ks = k(xs)

p2 = plt.figure()
plt.plot(ks,rs)
plt.xlim([0,40])
plt.ylim([0,0.8])

def model(x,r):
    a = x*(1+x**2)
    b = r**2*(1+x**2)-x
    return a/b
    
xs = np.linspace(0.1,20,n)
rs = np.linspace(0.1,1,n)
X,R = np.meshgrid(xs,rs)

K = model(X,R)

#fig3=plt.figure()
#ax2 = fig3.add_subplot(111, projection='3d')

#ax2.plot_surface(K, R, X)


 