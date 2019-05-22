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


n = 100000 # samples

# Timescales
alphax = random.uniform(1, 2, size=n) 
alphay = random.uniform(0, 1, size=n)
# Setup samples so that r is in [0,1]
r = alphay/alphax # turnover rate of predators (units multiples of prey)



# Beta parameters
d = random.uniform(size=n) # value of D*/X*
beta1 = d/alphax
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

# Condition for Saddle-Node, detJ = 0
saddle_sx = (a*b-c)/b

# Condition for Hopf, tr J = 0
hopf_sx = a -r*b



#fig1 = plt.figure()
#plt.scatter(range(n), s_x, color='red')
#plt.scatter(range(n), s_x2, color='blue')

# List of parameter names
lst = ['r', 'beta','Fx','Fy','Mx','My','Saddle Sx','Hopf Sx']

data = dict.fromkeys(lst,0)
data['r'] = r
data['beta'] = beta1
data['Fx'] = f_x
data['Fy'] = f_y
data['Mx'] = m_x
data['My'] = m_y
data['Saddle Sx'] = saddle_sx
data['Hopf Sx'] = hopf_sx



df = pd.DataFrame(data)

dfSaddle = df[abs(df['Saddle Sx']) <= 0.1 ]
dfHopf = df[abs(df['Hopf Sx']) <= 0.1 ]

#Plot pairs of parameters which give Saddle bifurcation
fig2 = plt.figure()
#plt.scatter(range(df2.shape[0]), df2['Sx'])
dfSaddle = dfSaddle.drop(columns='Saddle Sx')
g = sns.pairplot(dfSaddle, palette="husl")
g.fig.suptitle('Saddle Node Bifurcation Parameters')

#Plot pairs of parameters which give Hopf bifurcation
fig3 = plt.figure()
dfHopf = dfHopf.drop(columns='Hopf Sx')
g = sns.pairplot(dfSaddle, hue='species')
g.fig.suptitle('Hopf Bifurcation Parameters')
fig4 = plt.figure()
ax = fig4.add_subplot(111, projection='3d')
ax.scatter(dfHopf['beta'], dfHopf['r'], dfHopf['Fx'])
ax.set_xlabel('beta')
ax.set_ylabel('r')
ax.set_zlabel('fx')


