# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
from numpy import random
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# X = prey, Y = predators
# Timescales
alphax = 1/12. #months
alphay = 1. 

n = 10000 # samples
# Beta parameters
beta1 = random.uniform(size=n)
beta2 = 1 - beta1

# Elasticity parameters
s_x = random.uniform(0,3,size=n) # Prey birth elasticity
f_x = random.uniform(0,3,size=n) # Predation elasticity w resp. x
f_y = random.uniform(0,3,size=n) # Predation elasticity w resp y
m_x = random.uniform(0,3,size=n) # Mortality elasticiy of prey
m_y = random.uniform(0,3,size=n) # Mortality elasticiy of predators

# Smaller calculations
a = beta1*m_x+beta2*f_x
b = f_y-m_y
c = beta2*f_x*f_y

s_x2 = (a*b-c)/b

fig1 = plt.figure()
plt.scatter(range(n), s_x, color='red')
plt.scatter(range(n), s_x2, color='blue')

# List of parameter names
lst = ['beta','Fx','Fy','Mx','My','Sx']

data = dict.fromkeys(lst,0)
data['beta'] = beta1
data['Fx'] = f_x
data['Fy'] = f_y
data['Mx'] = m_x
data['My'] = m_y
data['Sx'] = s_x2

df = pd.DataFrame(data)

df2 = df[abs(df['Sx']) <= 0.1 ]
fig2 = plt.figure()
plt.scatter(range(df2.shape[0]), df2['Sx'])
df2 = df2.drop(columns='Sx')

g = sns.pairplot(df2, palette="husl")
