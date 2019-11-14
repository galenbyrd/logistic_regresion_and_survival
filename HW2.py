#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 10:32:47 2019

@author: GalenByrd
"""

import math
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.tools.tools as tools
import scipy.stats as stats

################# QUIESTION 1 ##############################################
df = pd.read_table('ICU/ICU.txt')
y = np.ravel(df[['STA']])

X = df[['INF']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())

table = np.array(pd.crosstab(df.STA,df.INF))
oddsratio, pvalue = stats.fisher_exact(table)
print("OddsR: ", oddsratio, "p-Value:", pvalue)
np.exp(results.params)


np.sqrt(np.sum(1/table))


np.exp(0.2074)
np.exp(1.6252)

################# QUIESTION 2 ##############################################
lowbwt = pd.read_table('LOWBWT/LOWBWT.txt')
y2 = np.ravel(lowbwt[['LOW']])

X2 = pd.get_dummies(lowbwt['RACE'],drop_first=True)
X2 = tools.add_constant(X2)
model = sm.Logit(y2,X2)
results = model.fit()
print(results.summary2())

table2 = np.array(pd.crosstab(lowbwt.RACE,lowbwt.LOW))
oddsratio, pvalue = stats.fisher_exact(table2[0:2])
print("OddsR: ", oddsratio, "p-Value:", pvalue)

oddsratio, pvalue = stats.fisher_exact(np.delete(table2,1, axis=0))
print("OddsR: ", oddsratio, "p-Value:", pvalue)

np.exp(results.params)



np.sqrt(np.sum(1/table2[0:2]))
np.sqrt(np.sum(1/np.delete(table2,1, axis=0)))



np.exp(-0.0635)
np.exp(1.7531)
np.exp(-0.0456)
np.exp(1.3179)


################# QUIESTION 3 ##############################################
df = pd.read_table('ICU/ICU.txt')
y = np.ravel(df[['STA']])

df['Interaction']=df['CRN']*df['AGE']

X = df[['CRN']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())
X = df[['CRN','AGE']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())
X = df[['CRN','AGE','Interaction']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())




################# QUIESTION 4 ##############################################
df = pd.read_table('BURN/BURN1000.txt')
y = np.ravel(df[['DEATH']])

df['Interaction']=df['INH_INJ']*df['AGE']

X = df[['INH_INJ']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())
X = df[['INH_INJ','AGE']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())
X = df[['INH_INJ','AGE','Interaction']]
X = tools.add_constant(X)
model = sm.Logit(y,X)
results = model.fit()
print(results.summary2())
covar=np.array(results.cov_params())

ages = [20,40,60,80]
oddsRatios = []
standErr = []
CI = []

for age in ages:
    oddsRatios.append(np.exp(6.2122+ (-0.0689*age)))
    standErr.append(np.sqrt(covar[1,1]+age**2*covar[3,3]+2*age*covar[1,3]))

for i,ratio in enumerate(oddsRatios):
    upper = np.log(ratio) + 1.96*standErr[i]
    lower = np.log(ratio) - 1.96*standErr[i]
    CI.append([np.exp(lower),np.exp(upper)])
    




###################### CODE GRAVE ########################################


df = pd.read_table('GLOW/GLOW500.txt')
y = np.ravel(df[['FRACTURE']])

df['Interaction']=df['PRIORFRAC']*df['AGE']
XA = df[['AGE','PRIORFRAC','Interaction']]
XA = tools.add_constant(XA)
model = sm.Logit(y,XA)
results = model.fit()
print(results.summary2())
covar=np.array(results.cov_params())

np.sqrt(covar[1,1]+55**2*covar[3,3]+2*55*covar[1,3])

np.exp(4.9613 + (-0.0574*55))
np.exp(4.9613 + (-0.0574*60))
np.exp(4.9613 + (-0.0574*65))
np.exp(4.9613 + (-0.0574*70))
np.exp(4.9613 + (-0.0574*75))
np.exp(4.9613 + (-0.0574*80))










