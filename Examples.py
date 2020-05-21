#!/usr/bin/env python
# coding: utf-8

# In[5]:


'''
Usage examples of option pricing for COS method, Filon's method and FFT method.
'''

from Ficamacos import model
from Ficamacos import method

## Case 1 ###

r = 0.05
T = 1
sigma = 0.1
S = 30
K = 20

param1 = model.GBMParam(sigma = sigma)
model1 = model.GBM(param1,r,T)
method.Filon(model1,S,K)
method.FFT(model1,S,K)





# In[ ]:


## Case 2 ###
r = 0.05
T = 1
sigma = 0.1
S = 40
K = 20
theta=0.05
v=0.1

param2 = model.VGParam(sigma = sigma, theta = theta, v = v)
model2 = model.VG(param2,r,T)
method.Filon(model2,S,K)
method.COS(model2,S,K,L=3)

