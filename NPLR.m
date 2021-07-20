classdef NPLR < matlab.mixin.Copyable
    %% Implements a Nested Peicewise Linear Representation (NPLR) of a function
    % obj = NPLR(N,Y,X) creates an NPLR object
    % obj = NPLR(N,Y,X,FLAG)
    % where
    %     N - No. of Dimensions
    %     Y - Data to be fitted, N-Dimensional Array. For vector-valued
    %         functions cell array of N-Dimensional arrays
    %     X - cell array of size N, i-th element contains vector x_i along
    %         which data is given.
    %  FLAG - [Optional] 'New' or 'Old' (switch between both implementations)
    %         default is 'New' (numerically much stable and fast especially
    %                           for vector-valued functions)
    %--------------------------------------------------------------------------
    properties (Access = public)
        N      % No of Dimensions
        nSize  % Size along each dimension
        nOut   % No of Outputs (Vector-Valued Functions)
        mu     % CPLR mu's
        Gamma  % NPLR Coeff Matrix (New Implementation)
        M      % CPLR Coeff Matrix (Old Implementation)
    end
    properties %(Access = private)
        Y      % Orignal Data Y
        X      % Orignal Data X
        Xhat   %
        
        OLD_FLAG % Implementation Flag
    end
    methods (Access = public)
        %-------------------------------------
        % Class Constructor
        function obj = NPLR(N,X,Y,varargin)
            validateattributes(N, {'numeric'},{'>',0});
            validateattributes(X, {'cell'},{'numel',N});
            validateattributes(Y, {'numeric','cell'},{});
            obj.OLD_FLAG = 0;
            if nargin > 3
                if strcmpi(varargin{1},'Old')
                    obj.OLD_FLAG = 1;
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
            % Compute mu's
            obj.mu = Compute_Mu(obj);
            % Compute Xhat
            obj.Xhat = cell(1,obj.N);
            for k = 1:obj.N
                n = obj.nSize(k);
                xvec = obj.X{k};
                obj.Xhat{k} = [ones(1,n);(xvec)';abs(ones(n-2,1)*(xvec)' - (obj.mu{k})*ones(n,1)')]';
            end
            
            % Compute Coeff. Matrix
            if obj.OLD_FLAG
                for i=1:obj.nOut
                    obj.M{i} = obj.CPLRFit_old(Y{i}); % Old Implementation
                end
            else
                obj.Gamma = obj.NPLRFit(Y); % New Implementation
            end
        end
        %-------------------------------------
        % N-D NPLR Implementation
        function Yq = eval(obj,Xq,varargin)
            % Evaluates NPLR or Compute Derivative (any first order or mixed
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
            
            if obj.OLD_FLAG
                nq  = size(Xq,1);
                Yq = zeros(nq,obj.nOut);
                for i=1:nq
                    Yq(i,:) = EvalOld(obj,Xq(i,:),DervFlag);
                end
            else
                Yq = EvalNew(obj,Xq,DervFlag);
            end
        end
        %-------------------------------------
        % Compute Multi-Linear Affine Approximation i.e Jacobian Matrix and Offset
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
            %             for k=1:obj.nOut
            %                 obj.Y{k} = obj1.Y{k} + obj2.Y{k};
            %             end
        end
    end
    methods (Access = private)
        %-------------------------------------
        % Compute Mu's
        function mu = Compute_Mu(obj)
            for i=1:obj.N
                mu{i} = obj.X{i}(2:end-1);
            end
        end
        %-------------------------------------
        % Fit ND NPLR (New Implementation)
        function [Gamma] = NPLRFit(obj,Y,ORDER)
            if nargin < 3
                ORDER = 1:obj.N;
            end
            I = flip(ORDER);
            Gamma = zeros(obj.nOut,prod(obj.nSize));
            for i=1:obj.nOut
                if obj.N > 1
                    GG = permute(Y{i},flip(ORDER));
                else
                    GG = Y{i};
                end
                for k = 1:obj.N
                    GG = reshape(GG,obj.nSize(I(k)),[]);
                    GG = ((obj.Xhat{I(k)})\GG)';
                end
                Gamma(i,:) = reshape(GG,1,[]);
            end
        end
        %-------------------------------------
        % Evaluate ND NPLR (New Implementation)
        function [Yq] = EvalNew(obj,Xq,DervFlag)
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
        % Fit ND CPLR (Old Implementation)
        function [M] = CPLRFit_old(obj,Y)
            % This Functions Computes CPLR Parameters for N-D Curve Fitting
            NSize = obj.nSize;
            NN = obj.N;
            if NN == 1
                % fit about x1
                M = (obj.Xhat{1})\Y;
            elseif NN == 2
                % fit about x1 and x2
                M = (obj.Xhat{1})\Y/(obj.Xhat{2})';
            elseif NN == 3
                % fit about x1 and x2
                for i = 1:NSize(3)
                    W2(:,:,i) = (obj.Xhat{1})\Y(:,:,i)/(obj.Xhat{2})';
                end
                % fit about x3
                Wk = W2;
                W2 = reshape(W2,[prod(NSize(1:NN-1)) NSize(NN)]);
                M = W2/(obj.Xhat{NN})';
            elseif NN >= 4
                % fit about x1 and x2
                ksize = NSize(3:end);
                OUT1 = reshape(Y,[NSize(1) NSize(2) prod(ksize)]);
                for i = 1:prod(ksize)
                    W2(:,:,i) = (obj.Xhat{1})\OUT1(:,:,i)/(obj.Xhat{2})';
                end
                Wk = W2;
                % fit about xk for 3 <= k <= N-1
                for k=3:NN-1
                    Wk = reshape(Wk,[prod(NSize(1:k-1)) NSize(k) prod(NSize(k+1:NN))]);
                    temp = zeros(size(Wk));
                    for i = 1:prod(NSize(k+1:NN))
                        temp(:,:,i) = Wk(:,:,i)/(obj.Xhat{k})';
                    end
                    Wk = temp;
                end
                % fit about xN
                Wk = reshape(Wk,[prod(NSize(1:NN-1)) NSize(NN)]);
                M  = Wk/(obj.Xhat{NN})';
            end
        end
        %-------------------------------------
        % Evaluate ND CPLR (Old Implementation)
        function [Yq] = EvalOld(obj,xq,DervFlag)
            Yq = zeros(1,obj.nOut);
            for iO=1:obj.nOut
                W = obj.M{iO};
                if obj.N == 1
                    x1q = xq(1);
                    x1vec = obj.HAT(x1q,obj.mu{1},DervFlag(1));
                    Yq = x1vec'*W;
                elseif obj.N >= 2
                    % 1D CPLR w.r.t xq_k for 3 <= k <= N
                    for k = obj.N:-1:3
                        xkq   = xq(k);
                        xkvec = obj.HAT(xkq,obj.mu{k},DervFlag(k));
                        ksize = [prod(obj.nSize(1:k-2)) obj.nSize(k-1)];
                        W     = reshape(W*xkvec,ksize);
                    end
                    % 2D CPLR w.r.t xq_1 & xq_2
                    x1q = xq(1); x2q = xq(2);
                    x1vec = obj.HAT(x1q,obj.mu{1},DervFlag(1));
                    x2vec = obj.HAT(x2q,obj.mu{2},DervFlag(2));
                    Yq(iO) = x1vec'*W*x2vec;
                end
            end
        end
        %-------------------------------------
    end
    methods (Static)
        %-------------------------------------
        % Compute OrderReversing Indexes
        function [Ir] = IndexReversedOrder(nSize)
            Sv = @(p,q)(reshape(reshape(uint32(1:p*q),q,[])',[],1));
            Ir = Sv(prod(nSize(2:end)),nSize(1));
            for k=2:length(nSize)-1
                m  = prod(nSize(1:k-1));
                t_ = Sv(prod(nSize(k+1:end)),nSize(k));
                pq = prod(nSize(k+1:end))*nSize(k);
                tempv = kron(uint32(ones(m,1)),t_) + uint32(kron((0:m-1)',pq*ones(pq,1)));
                Ir = Ir(tempv);
            end
        end
        %-------------------------------------
        % Compute Shuffling Matrix
        function [S] = ShufflingMatrix(p,q)
            r = p*q;
            I = eye(r);
            S = I(1:q:r,:);
            for k=2:q
                S = [S;I(k:q:r,:)];
            end
        end
        %-------------------------------------
    end
    methods (Access = private, Static)
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
    end
end