"""
Created on Fri Apr  9 21:53:12 2021

@author: Zeeshan Iqbal Khan
"""
from numpy.random import random as rand
import numpy as np
from NPLRTools import NPLR, NPLR_TEST_DATA



N  = 2  # Number of Dimensions
nO = 2  # Number of OUtputs

X,Y,nSize = NPLR_TEST_DATA(N,nO)

obj = NPLR(N,X,Y)

objN = NPLR(Gamma=obj.Gamma,mu=obj.mu)

# %% Generate Refined Test Data (within the limits)
nq = 10;
Xq = []
for k in np.arange(0,N):
    Xq.append(rand([nq,1])*(max(X[k]) - min(X[k])) + min(X[k]))

Xq = np.hstack(Xq)

Yq = obj.Eval(Xq)

B,eps,Yq_ = obj.Jacobian(Xq)

mdict = {'Xq':Xq, 'Yq':Yq, 'B':B, 'eps':eps, 'Yq_':Yq_}

#%% Save data for validation in matlab
obj.SaveData('temp.mat',mdict)

