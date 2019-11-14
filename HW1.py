#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 09:46:03 2019

@author: GalenByrd
"""
import math
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.tools.tools as tools

# IMPORT DATA
df = pd.read_table('GLOW/GLOW500.txt')
y = np.ravel(df[['FRACTURE']])

# 1A
XA = df[['AGE','WEIGHT','PRIORFRAC','PREMENO']]
encoding = pd.get_dummies(df['RATERISK'],drop_first=True)
XA=pd.concat([XA,encoding], axis=1)
XA = tools.add_constant(XA)

model = sm.Logit(y,XA).fit()
print(model.summary2())

XB = df[['AGE','PRIORFRAC']]
XB=pd.concat([XB,encoding], axis=1)
XB = tools.add_constant(XB)

model2 = sm.Logit(y,XB).fit()
print(model2.summary2())

# 1 B
#G = -2[Log-Likelihood - Log-Likelihood]
#G = -2[-259.45-(-259.04)]=0.8324

# 2 D
df = pd.read_table('ICU/ICU.txt')
y = np.ravel(df[['STA']])

XA = df[['AGE','CAN','CPR','INF']]
encoding = pd.get_dummies(df['RACE'],drop_first=True)
XA=pd.concat([XA,encoding], axis=1)
XA = tools.add_constant(XA)
model = sm.Logit(y,XA)
results = model.fit()
print(results.summary2())


XB = df[['AGE','CPR']]
XB = tools.add_constant(XB)
model2 = sm.Logit(y,XB)
results2 = model2.fit()
print(results2.summary2())

cov = df[['AGE','CPR','STA']].vcov()
results2.cov_params()

const = 200*[1]
logit_const = sm.Logit(df['STA'],const)

result_const = logit_const.fit().summary2()
print(result_const)


ob = np.matrix([1,60,1])
math.sqrt(ob*np.matrix(results2.cov_params())*np.transpose(ob))


