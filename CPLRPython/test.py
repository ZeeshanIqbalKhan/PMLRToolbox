# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 22:05:34 2020

@author: Zeeshan Iqbal Khan
"""
import numpy as np
import cplr

N = 1;

nSize = np.random.randint(1,10,N)

OUT = np.random.random_sample(nSize)

obj = cplr.CPLR(1,OUT,X)

for k in range(N):
    n = nSize[k]
#     dx   = round(abs(rand*10 + 1));
#     x0   = round(abs(rand*10 + 1));
#     xf   = x0+dx*(n-1);
#     xvec = (x0:dx:xf)';
#     X{k} = xvec;