#!/usr/bin/env python
# coding: utf-8



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
print(method.Filon(model1,S,K))
print(method.FFT(model1,S,K))



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
print(method.Filon(model2,S,K))
print(method.COS(model2,S,K,L=3))

