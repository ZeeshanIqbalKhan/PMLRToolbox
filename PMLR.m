classdef PMLR < matlab.mixin.Copyable
    %% Implements a Peicewise Multi-Linear Representation (PMLR) of a function
    % obj = PMLR(N,Y,X) creates a PMLR object
    % obj = PMLR(N,Y,X,FLAG)
    % where
    %     N - No. of Dimensions
    %     Y - Data to be fitted, N-Dimensional Array. For vector-valued
    %         functions cell array of N-Dimensional arrays
    %     X - cell array of size N, i-th element contains vector x_i along
    %         which data is given.
    %  FLAG - [Optional] 'BKRON' or 'OTHER'(switch between both methods to compute Gamma)
    %         default is 'BKRON' as derived in papers  
    %                    'OTHER' is (numerically faster especially for vector-valued functions)
    %                                
    %--------------------------------------------------------------------------
    properties (Access = public)
        N      % No of Dimensions
        nSize  % Size along each dimension
        nOut   % No of Outputs (Vector-Valued Functions)
        mu     % PMLR mu's (interior grid points)
        Gamma  % PMLR Coeff Matrix
    end
    properties (Access = private)
        Y       % Orignal Data Y
        X       % Orignal Data X
        Xhat    % square matrices of xhats at each grid point
        XhatInv % Inverse of Xhat
        GAMMA_FLAG % PMLR Gamma Computation Algo [1 for BKRON, 0 for OTHER]
    end
    methods (Access = public)
        %-------------------------------------
        % Class Constructor
        function obj = PMLR(N,X,Y,varargin)
            validateattributes(N, {'numeric'},{'>',0});
            validateattributes(X, {'cell'},{'numel',N});
            validateattributes(Y, {'numeric','cell'},{});
            obj.GAMMA_FLAG = 1;
            if nargin > 3
                if strcmpi(varargin{1},'BKRON')
                    obj.GAMMA_FLAG = 1;
                elseif strcmpi(varargin{1},'OTHER')
                    obj.GAMMA_FLAG = 0;
                end
            end
            
            obj.N = uint8(N);
            
            if iscell(Y)
                obj.nOut = uint8(numel(Y));
                temp = cell2mat(cellfun(@(x)size(x)',Y,'UniformOutput',false));
                if norm(temp - repmat(temp(:,1),1,size(temp,2))) ~= 0
                    error('For vector-valued functions, size of each ND array must be same.')
                end
                SS = cellfun(@(z)length(z),X);
                if numel(SS)==1, SS = [SS 1]; end
                validateattributes(Y{1}, {'numeric'},{'size',SS});
                obj.nSize = size(Y{1});
            else
                validateattributes(Y, {'numeric'},{'size',cellfun(@(z)length(z),X)});
                obj.nOut = uint8(1);
                obj.nSize = size(Y);
            end
            obj.Y = Y;
            obj.X = X;
            % Extract mu's
            obj.mu = getMu(obj);
            % Compute Xhat
            obj.Xhat = cell(1,obj.N);
            obj.XhatInv = cell(1,obj.N);
            for k = 1:obj.N
                n = obj.nSize(k);
                xvec = obj.X{k};
                obj.Xhat{k} = [ones(1,n);(xvec)';abs(ones(n-2,1)*(xvec)' - (obj.mu{k})*ones(n,1)')];
                obj.XhatInv{k} = eye(n)/obj.Xhat{k};
            end
            
            % Compute Coeff. Matrix
            if obj.GAMMA_FLAG
                obj.Gamma = obj.PMLR_BKronFit(Y); % BKRON Implementation
            else
                obj.Gamma = obj.PMLRFit(Y); % Other Implementation
            end
        end
        %-------------------------------------
        % N-D PMLR Implementation
        function Yq = eval(obj,Xq,varargin)
            % Evaluates PMLR or Compute Derivative (any first order or mixed
            % multi-order derivaties i.e. all derivatives which exits for
            % piece-wise multi-linear functions
            % e.g. d^3/dx1dx2dx3 exits while d^2/dx1^2 doesnt
            p = inputParser;
            addRequired(p,'Xq',@(x)(size(x,2) == obj.N));
            addParameter(p,'DerivativeIndices',[],@(x)(isnumeric(x) && length(unique(x)) <= obj.N));
            parse(p,Xq,varargin{:});
            
            Xq = p.Results.Xq;
            D_ind = unique(p.Results.DerivativeIndices);
            DervFlag = zeros(1,obj.N);
            for k=1:obj.N
                DervFlag(k) = ismember(k,D_ind);
            end
            
            Yq = obj.EvalPMLR(Xq,DervFlag);
        end
        %-------------------------------------
        % Compute Local Affine Approximation i.e Jacobian Matrix and Offset
        function [B,eps,Yq] = Jacobian(obj,Xq,varargin)
            p = inputParser;
            addRequired(p,'Xq',@(x)(size(x,2) == obj.N));
            addParameter(p,'DerivativeIndices',num2cell(1:obj.N),...
                @(x)(all(cellfun(@(y)(isnumeric(y) && length(unique(y)) <= obj.N),x))));
            parse(p,Xq,varargin{:});
            
            Xq = p.Results.Xq;
            D_ind = cellfun(@(x)unique(x),p.Results.DerivativeIndices,'UniformOutput',false);
            
            m = length(D_ind);
            Yq = reshape(obj.eval(Xq)',obj.nOut,1,[]);
            B = zeros(obj.nOut,m,size(Xq,1));
            eps = Yq;
            for i = 1:m
                B(:,i,:) = obj.eval(Xq,'DerivativeIndices',D_ind{i})';
                XX = reshape(repmat(prod(Xq(:,D_ind{i}),2)',obj.nOut,1),obj.nOut,1,[]);
                eps = eps - B(:,i,:).*XX;
            end
        end
        %-------------------------------------
        function obj = plus(obj1,obj2)
            assert(obj1.N == obj2.N);
            assert(obj1.nOut == obj2.nOut);
            assert(all(obj1.nSize == obj2.nSize));
            for k=1:obj1.N
                assert(all(obj1.mu{k} == obj2.mu{k}));
            end
            obj = obj1.copy();
            obj.Gamma = obj1.Gamma + obj2.Gamma;
            for k=1:obj.nOut
                obj.Y{k} = obj1.Y{k} + obj2.Y{k};
            end
        end
        %-------------------------------------        
    end
    methods (Access = private)
        %-------------------------------------
        % Extract Mu's (interior grid points)
        function mu = getMu(obj)
            for i=1:obj.N
                mu{i} = obj.X{i}(2:end-1);
            end
        end
        %-------------------------------------
        % Fit ND PMLR (BKron Implementation)
        function [Gamma] = PMLR_BKronFit(obj,Y,I)
            if nargin < 3
                I = 1:obj.N;
            end
            Gamma = zeros(obj.nOut,prod(obj.nSize));
            for i=1:obj.nOut
                if obj.N > 1
                    GG = permute(Y{i},flip(I));
                else
                    GG = Y{i};
                end
                GG = vec(GG);
                
                for k=1:obj.N
                    kk_ = prod(obj.nSize(I(k+1:end)));
                    GG = reshape(GG,kk_,[]);
                    GG = obj.bkron(GG,obj.XhatInv{I(k)});
                end
                Gamma(i,:) = GG;
            end
        end
        %-------------------------------------
        % Fit ND PMLR (Other Implementation)
        function [Gamma] = PMLRFit(obj,Y,ORDER)
            if nargin < 3
                ORDER = 1:obj.N;
            end
            I = flip(ORDER);
            Gamma = zeros(obj.nOut,prod(obj.nSize));
            for i=1:obj.nOut
                if obj.N > 1
                    GG = permute(Y{i},I);
                else
                    GG = Y{i};
                end
                for k = 1:obj.N
                    GG = reshape(GG,obj.nSize(I(k)),[]);
                    GG = (obj.XhatInv{I(k)}'*GG)';
                end
                Gamma(i,:) = reshape(GG,1,[]);
            end
        end
        %-------------------------------------
        % Evaluate ND PMLR 
        function [Yq] = EvalPMLR(obj,Xq,DervFlag)
            nq  = size(Xq,1);
            Yq = zeros(nq,obj.nOut);
            for i=1:nq
                x1q = Xq(i,1);
                x1vec = obj.HAT(x1q,obj.mu{1},DervFlag(1));
                temp = x1vec;
                for k = 2:obj.N
                    xkq = Xq(i,k);
                    xkvec = obj.HAT(xkq,obj.mu{k},DervFlag(k));
                    temp = kron(temp,xkvec);
                end
                Yq(i,:) = obj.Gamma*temp;
            end
        end
        %-------------------------------------
    end
    methods (Static, Access = private) 
        %-------------------------------------
        % Hat Operator
        function xhat = HAT(x,mu,isDerivative)
            if(~isDerivative)
                xhat = [1;x;abs(x - mu)];
            else
                xhat = [0;1;sign(x - mu)]; % d/dx(xhat)
            end
        end
        %-------------------------------------
        % Block Kronecker
        function [C] = bkron(A,B)
            [m,n] = size(A);
            [p,q] = size(B);
            % assert(mod(n,p) == 0);
            k = n/p;
            if(1)
                Ah = reshape(A',p,[])';
                II = reshape(reshape(1:m*k,k,[])',1,[]);
                C = reshape(Ah(II,:)*B,m,[]);
                II = reshape(reshape(1:q*k,k,[])',1,[]);
                C = C(:,II);
            else
                C = A*kron(speye(k),B);
            end
        end
        %-------------------------------------
    end
end