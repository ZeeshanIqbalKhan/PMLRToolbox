classdef CPLR
    %% Implements a Canonical Peicewise Linear Representation (CPLR) of a function
    % obj = CPLR(N,Y,X) creates an CPLR object, where
    %     N - No. of Dimensions
    %     Y - Data to be fitted, N-Dimensional Array
    %     X - cell array of size N, i-th element contains vector x_i along
    %         which data is given.
    %--------------------------------------------------------------------------
    properties (Access = public)
        N      % No of Dimensions
        nSize  % Size along each dimension
        M      % CPLR Coeff Matrix
        mu     % CPLR mu's
    end
    properties (Access = private)
        Y      % Orignal Data Y
        X      % Orignal Data X
    end
    methods (Access = public)
        %-------------------------------------
        % Class Constructor
        function obj = CPLR(N,Y,X)
            obj.N = N;
            obj.nSize = size(Y);
            obj.Y = Y;
            obj.X = X;
            obj.mu = Compute_Mu(obj);
            for k = 1:N
                n = obj.nSize(k);
                xvec = obj.X{k};
                Xhat{k} = [ones(1,n);(xvec)';abs(ones(n-2,1)*(xvec)' - (obj.mu{k})*ones(n,1)')]';
            end
            obj.M = CPLRFit(obj,Xhat);
        end
        %-------------------------------------
        % Evaluate CPLR
        function [Yq,dYq,eps] = eval(OBJ,Xq)
            % Evaluate CPLR and/or compute slopes and intercepts
            % [Yq,dYq,eps] = eval(OBJ,Xq) computes ouputs of CPLR at given query points
            % Inputs:
            %    OBJ: CPLR object or 1D array of CPLR objects (each object
            %         must have same grid points and dimensions)
            %     Xq: N x nq matrix of query points, where N must be same
            %         as dimension of OBJ (OBJ.N)
            %
            % Outputs:
            %     Yq: Nobj x nq matrix of output of CPLR, (where Nobj is 
            %         number of CPLR objects in OBJ)
            %    dYq: [Optional] Nobj x N x nq matrix, Slopes matrix at
            %                    each query point
            %    eps: [Optional] Nobj x nq matrix, Intercept vector at each
            %                    query point
            %
            Nobj = length(OBJ);
            if(~all(size(Xq,1) == [OBJ.N]))
                error('Size Mismatch');
            end
            nq  = size(Xq,2);
            Yq  = zeros(Nobj,1,nq);
            %--------------------------------------------------------------
            % Compute CPLR
            for KK = 1:Nobj
                obj = OBJ(KK);
                % N-D CPLR Implementation
                for i=1:nq
                    if obj.N == 1
                        x1q = Xq(1,i);
                        x1vec = [1;x1q;abs(x1q - obj.mu{1})];
                        Yq(KK,1,i) = x1vec'*obj.M;
                    elseif obj.N >= 2
                        % 1D CPLR w.r.t xq_k for 3 <= k <= N
                        W = obj.M;
                        for k = obj.N:-1:3
                            xkq   = Xq(k,i);
                            xkvec = [1;xkq;abs(xkq - obj.mu{k})];
                            ksize = [prod(obj.nSize(1:k-2)) obj.nSize(k-1)];
                            W     = reshape(W*xkvec,ksize);
                        end
                        % 2D CPLR w.r.t xq_1 & xq_2
                        x1q = Xq(1,i); x2q = Xq(2,i);
                        x1vec = [1;x1q;abs(x1q - obj.mu{1})];
                        x2vec = [1;x2q;abs(x2q - obj.mu{2})];
                        Yq(KK,1,i) = x1vec'*W*x2vec;
                    end
                end
            end
            
            %--------------------------------------------------------------
            % Compute Slopes & Intercept
            if nargout > 1
                %----------------------------------------------------------
                % Compute Slopes
                nX  = size(Xq,1);
                dYq = zeros(Nobj,nX,nq);
                for JJ = 1:nX
                for KK = 1:Nobj
                    obj = OBJ(KK);
                    % N-D CPLR Implementation
                    for i=1:nq
                        if obj.N == 1
                            x1q = Xq(1,i);
                            if JJ==1
                                x1vec = [0;1;sign(x1q - obj.mu{1})];
                            else
                                x1vec = [1;x1q;abs(x1q - obj.mu{1})];
                            end
                            dYq(KK,1,i) = x1vec'*obj.M;
                        elseif obj.N >= 2
                            % 1D CPLR w.r.t xq_k for 3 <= k <= N
                            W = obj.M;
                            for k = obj.N:-1:3
                                xkq   = Xq(k,i);
                                if JJ == k
                                    xkvec = [0;1;sign(xkq - obj.mu{k})];
                                else
                                    xkvec = [1;xkq;abs(xkq - obj.mu{k})];
                                end
                                ksize = [prod(obj.nSize(1:k-2)) obj.nSize(k-1)];
                                W     = reshape(W*xkvec,ksize);
                            end
                            % 2D CPLR w.r.t xq_1 & xq_2
                            x1q = Xq(1,i); x2q = Xq(2,i);
                            if JJ == 1
                                x1vec = [0;1;sign(x1q - obj.mu{1})];
                            else
                                x1vec = [1;x1q;abs(x1q - obj.mu{1})];
                            end
                            if JJ == 2
                                x2vec = [0;1;sign(x2q - obj.mu{2})];
                            else
                                x2vec = [1;x2q;abs(x2q - obj.mu{2})];
                            end
                            dYq(KK,JJ,i) = x1vec'*W*x2vec;
                        end
                    end
                end
                end
                %----------------------------------------------------------
                % Compute Intercepts
                eps = zeros(Nobj,1,nq);
                for i = 1:nq
                    eps(:,:,i) = Yq(:,:,i) - dYq(:,:,i)*Xq(:,i);
                end
                eps = squeeze(eps);
            end
            Yq = squeeze(Yq);
        end
        %-------------------------------------
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
        % Fit ND CPLR
        function [M] = CPLRFit(obj,Xhat)
            % This Functions Computes CPLR Parameters for N-D Curve Fitting
            NSize = obj.nSize;
            NN = obj.N;
            if NN == 1
                % fit about x1
                M = (Xhat{1})\obj.Y;
            elseif NN == 2
                % fit about x1 and x2
                M = (Xhat{1})\obj.Y/(Xhat{2})';
            elseif NN == 3
                % fit about x1 and x2
                for i = 1:NSize(3)
                    W2(:,:,i) = (Xhat{1})\obj.Y(:,:,i)/(Xhat{2})';
                end
                % fit about x3
                Wk = W2;
                W2 = reshape(W2,[prod(NSize(1:NN-1)) NSize(NN)]);
                M = W2/(Xhat{NN})';
            elseif NN >= 4
                % fit about x1 and x2
                ksize = NSize(3:end);
                OUT1 = reshape(obj.Y,[NSize(1) NSize(2) prod(ksize)]);
                for i = 1:prod(ksize)
                    W2(:,:,i) = (Xhat{1})\OUT1(:,:,i)/(Xhat{2})';
                end
                Wk = W2;
                % fit about xk for 3 <= k <= N-1
                for k=3:NN-1
                    Wk = reshape(Wk,[prod(NSize(1:k-1)) NSize(k) prod(NSize(k+1:NN))]);
                    temp = zeros(size(Wk));
                    for i = 1:prod(NSize(k+1:NN))
                        temp(:,:,i) = Wk(:,:,i)/(Xhat{k})';
                    end
                    Wk = temp;
                end
                % fit about xN
                Wk = reshape(Wk,[prod(NSize(1:NN-1)) NSize(NN)]);
                M  = Wk/(Xhat{NN})';
            end
        end
        %-------------------------------------
    end
end