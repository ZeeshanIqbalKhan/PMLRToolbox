clearvars; %clc;
rng(10) % For Repeatability
F = @(x)(exp(0.1*x) - x.*x.*x + 0.2*sin(x));
dF = @(x)(0.1*exp(0.1*x) - 3.*x.*x + 0.2*cos(x));
x0   = 0;    xf   = 10;
%% Data for Fitting & Testing
nO = 1;
nX = 1000;
N = 1; % Currently for only N=1, this file needs to be updated for N > 1

tic
nSize = nX*ones(N,1);
X = cell(1,N);
for k = 1:N
    n    = nSize(k);
    dx   = (xf-x0)/(n-1); % Make sure equal spacing 
    X{k} = round((x0:dx:xf)',10); % xvec
end
clear k n dx
Y = cell(1,nO);
for iO = 1:nO
        Y{iO} = F(X{1});
end
%% Generate Refined Test Data (within the limits)
nq = nX;
Xq = zeros(nq,N);
for k = 1:N
    Xq(:,k) = round(linspace(0,10,nq),10);
end
%% MATLAB Linear Interpolation
for iO=1:nO
    dx = 0.0001;
    dYinterp(:,iO) = (interp1(X{1},Y{1},min(Xq+dx,xf))-interp1(X{1},Y{1},max(Xq-dx,x0)))/(2*dx);
end
%% NPLR Model
obj = PMLR(N,X,Y);   % Create Object
dYpmlr  = obj.Jacobian(Xq); % Evaluate Object
%%
RMSE = sqrt(mean((dYinterp(:)-dYpmlr(:)).^2))/sqrt(mean(dYinterp(:).^2));
fprintf('RMSE = %10.4e for N = %2g in %7.4f (sec)\n',RMSE,N,toc);
%%
figure(1),clf
df = dF(Xq);
plot(Xq(2:end-1),squeeze(dYpmlr(2:end-1))-df(2:end-1),'.-b',...
     Xq(2:end-1),squeeze(dYinterp(2:end-1))-df(2:end-1),'.-r'),grid on
legend('PMLR vs. Analytical','Interpolation vs. Analytical')
title('Error in Jacobian')