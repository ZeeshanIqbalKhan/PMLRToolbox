# Implements a Canonical Peicewise Linear Representation (CPLR) of a function
# obj = CPLR(N,Y,X) creates an CPLR object, where
#     N - No. of Dimensions
#     Y - Data to be fitted, N-Dimensional Array
#     X - cell array of size N, i-th element contains vector x_i along
#         which data is given.
# 
# by Hafiz Zeeshan Iqbal Khan
import numpy as np
import matplotlib.pyplot as plt
class CPLR:
    # N      % No of Dimensions
    #     nSize  % Size along each dimension
    #     M      % CPLR Coeff Matrix
    #     mu     % CPLR mu's
    # end
    # properties (Access = private)
    #     Y      % Orignal Data Y
    #     X      % Orignal Data X
    # end
    # methods (Access = public)
    def __init__(self, N, Y, X):
        self.N = N;
        self.nSize = np.shape(Y);
        self.Y = Y;
        self.X = X;
        #self.mu = self.Compute_Mu(obj);
        #for k in range(N)
        #    n = obj.nSize(k);
        #    xvec = obj.X{k};
        #    Xhat{k} = [ones(1,n);(xvec)';abs(ones(n-2,1)*(xvec)' - (obj.mu{k})*ones(n,1)')]';
        #end
        #obj.M = CPLRFit(obj,Xhat);
        
        
        
X = [np.random.uniform(-1,1,10)]
Y = [np.random.uniform(-1,1,10)]   
plt.plot(X,Y)     
obj = CPLR(1,Y,X)