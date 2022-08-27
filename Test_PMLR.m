clearvars; %clc;
rng(10) % For Repeatability
%% Data for Fitting & Testing
nq = 100; % Number of query data points
nO = 3;    % Number of outputs
tic
for N = 1:10
nSize = round(abs(rand(1,N))*2 + 3); % Random # of grid points along each dim
Y = cell(1,nO);
for iO = 1:nO
    if N > 1
        Y{iO} = rand(nSize);
    else
        Y{iO} = rand(nSize,1);
    end
end
X = cell(1,N);
for k = 1:N
    n    = nSize(k);
    dx   = abs(rand*10 + 1); % Make sure equal spacing 
    x0   = abs(rand*10 + 1);
    xf   = x0+dx*(n-1);
    X{k} = (x0:dx:xf)'; % xvec
end
clear k n x0 dx xf
%% Generate Refined Test Data (within the limits)
Xq = zeros(nq,N);
for k = 1:N
    Xq(:,k) = rand(nq,1).*(max(X{k}) - min(X{k})) + min(X{k});
end
%% MATLAB Linear Interpolation
XQ = mat2cell(Xq,nq,ones(1,N));
[GRID{1:numel(X)}] = ndgrid(X{:});
for iO=1:nO
Yinterp(:,iO) = interpn(GRID{:},Y{iO},XQ{:});
end
%% PMLR Model
obj(N) = PMLR(N,X,Y,'BKRON'); % Create Object
Ypmlr  = obj(N).eval(Xq);   % Evaluate Object
%% Display Results
RMSE(N) = sqrt(mean((Yinterp(:)-Ypmlr(:)).^2))/sqrt(mean(Yinterp(:).^2));
fprintf('RMSE = %10.4e for N = %2g in %7.4f (sec)\n',RMSE(N),N,toc);
tic
end

%% Plot
figure(2)
semilogy(1:N,RMSE,'.-b','MarkerSize',14),xlabel('Dimensions - N'), ylabel('RMSE')
grid on,box on