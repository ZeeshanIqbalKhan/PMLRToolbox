"""
Created on Fri Apr  9 21:43:28 2021

@author: Zeeshan Iqbal Khan
"""
import numpy as np
from numpy.matlib import repmat
from numpy.random import random as rand
from scipy.io import savemat

def merge_dicts(*dict_args):
    """
    Given any number of dictionaries, shallow copy and merge into a new dict,
    precedence goes to key-value pairs in latter dictionaries.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


#%% Function to generate random data for NPLR testing
def NPLR_TEST_DATA(N,nO,Size=np.array([]),seed=-1234):
    if(seed != -1234):
        np.random.seed(seed) # seed for repeatability
        
    # Ensure interger sizes
    if(len(Size)==0):
        Size = np.abs(rand(N))*5 + 3

    nSize = []
    for a in Size:
        if(a<3):
            print("WARNING: size along any dimension cannot be less than 3, so converting it to 3")
        nSize.append(max(int(a),3))
    
    # Generate Y Random Data
    Y = []
    for iO in np.arange(0,nO):
        if(N==1):
            Y.append(np.reshape(rand(nSize),(-1,1)))
        else:
            Y.append(rand(nSize))
        
    X = []    
    for k in np.arange(0,N):
        dx   = abs(rand()*5 + 1) # Make sure equal spacing 
        x0   = rand()*10 + 1
        xf   = x0 + dx*(nSize[k]-1)
        X.append(np.reshape(np.arange(x0,xf,dx-0.001),(-1,1)))
    
    return X,Y,nSize

#%% Nested Piecewise Linear Representation
class NPLR():
    
    def __init__(self, *args,**kwargs): #N,X,Y)
        if(len(args)==3):
            N,X,Y = args
            self.DATA_FLAG = True
            self.N = N
            self.X = X
            self.Y = Y
            self.nOut = len(self.Y)
            self.nSize = np.shape(self.Y[0])
            self.mu, self.Xhat = self.Compute_MuandXhat()
            self.Gamma = self.NPLRFit()
        elif('Gamma' and 'mu' in kwargs.keys()):
            self.DATA_FLAG = False
            self.Gamma = kwargs.get('Gamma')
            self.mu = kwargs.get('mu')
            self.N  =  len(self.mu)
            self.nOut = self.Gamma.shape[0]
            self.nSize = []
            for a in self.mu:
                self.nSize.append(len(a))
            self.X = None
            self.Y = None
            self.Xhat = None
        
        

    
    def Compute_MuandXhat(self):
        if(not self.DATA_FLAG):
            ValueError("Cannot Call 'Compute_MuandXhat' because NPLR object dosent have data.")
        Xhat = []
        mu = []
        for i in np.arange(0,self.N):
            n = self.nSize[i]
            mu_ = self.X[i][1:-1]
            a_ = np.ones([1,n])
            b_ = np.transpose(self.X[i])
            c_ = abs(np.ones([n-2,1])@(np.transpose(self.X[i])) - (mu_)@np.ones([1,n]))
            mu.append(mu_)
            Xhat.append(np.transpose(np.vstack((a_,b_,c_))))
        return mu, Xhat
    
    def NPLRFit(self):
        if(not self.DATA_FLAG):
            ValueError("Cannot Call 'NPLRFit' because NPLR object dosent have data.")
        Gamma = [] #Gamma = np.zeros([self.nOut,np.prod(self.nSize)])
        for i in np.arange(0,self.nOut):
            GG = np.reshape(self.Y[i],(self.nSize[0],-1),order="F")
            GG = np.transpose(np.linalg.inv(self.Xhat[0])@GG)
            for k in np.arange(1,self.N):
                GG = np.reshape(GG,(self.nSize[k],-1),order="F")
                GG = np.transpose(np.linalg.inv(self.Xhat[k])@GG)
            Gamma.append(np.reshape(GG,(1,-1),order="F"))
            
        return np.vstack(Gamma)
    
    def Eval(self, Xq, D_ind=np.array([])):
        DervFlag = []
        for k in np.arange(0,self.N):
            DervFlag.append(np.any(np.isin(D_ind,k)))
        DervFlag = np.hstack(DervFlag)
        
                
        nq = Xq.shape[0]
        Yq = []
        for i in np.arange(0,nq):
            x1q = Xq[i,0]
            x1vec = self.HAT(x1q,self.mu[0],DervFlag[0])
            temp = x1vec
            for k in np.arange(1,self.N):
                xkq = Xq[i,k]
                xkvec = self.HAT(xkq,self.mu[k],DervFlag[k])
                temp = np.kron(xkvec,temp)   
            Yq.append(np.transpose(self.Gamma@temp))
        return np.vstack(Yq)            
    
    def Jacobian(self,Xq,DerivativeIndices=None):
        if(DerivativeIndices==None):
            DerivativeIndices = np.arange(0,self.N).tolist()
        
        _Dind = []
        for x in DerivativeIndices:
            _Dind.append(np.unique(x))
        
        m = len(_Dind)
        Yq = np.reshape(np.transpose(self.Eval(Xq)),(self.nOut,1,-1),order="F")
        B = np.zeros((self.nOut,m,Xq.shape[0]))
        eps = Yq
        for i in np.arange(0,m):
            B[:,i,:] = np.transpose(self.Eval(Xq,D_ind=_Dind[i]))
            _X = np.reshape(repmat(np.transpose(np.prod(Xq[:,_Dind[i]],axis=1)),self.nOut,1),(self.nOut,1,-1),order="F")
            eps = eps - np.reshape(B[:,i,:],(self.nOut,1,-1),order="F")*_X
        
        return B, eps, Yq

    
    def HAT(self,x,mu,isDerivative=0):
        if(isDerivative):            
            return np.vstack(([0],[1],np.sign(x - mu))) # d/dx(xhat)
        else:
            return np.vstack(([1],[x],abs(x - mu)))  # xhat
    
    def SaveData(self,path,m_dict):
        if(self.DATA_FLAG):
            _X = dict()
            for k in np.arange(0,self.N):
                _X['x_' + str(k)] = self.X[k]
    
            _Y = dict()
            for k in np.arange(0,self.nOut):
                _Y['y_' + str(k)] = self.Y[k]

            mdict = {'N': self.N, 'X': _X, 'Y': _Y, 'Gamma': self.Gamma}
        else:
            mdict = {'Gamma': self.Gamma}
            
        mdictn = merge_dicts(mdict, m_dict)
        savemat(path, mdictn)
    