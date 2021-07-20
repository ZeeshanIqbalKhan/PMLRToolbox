%clc; 
clear;
clearvars;
rng(0)
%% Data for Fitting & Testing
figure(1),clf
tic
nOut = 3;
for N = 1:10;
nSize = round(abs(rand(1,N))*2 + 3);
clear OUT
for kk=1:nOut
    if N > 1
        OUT{kk} = rand(nSize);
    else
        OUT{kk} = rand(nSize,1);
    end
end
for k = 1:N
    n    = nSize(k);
    dx   = round(abs(rand*10 + 1));
    x0   = round(abs(rand*10 + 1));
    xf   = x0+dx*(n-1);
    xvec = (x0:dx:xf)';
    X{k} = xvec;
end
clear k n x0 dx xf xvec mu
%% Generate Refined Test Data (within the limits)
nq = 1000;
for k = 1:N
    Xq(:,k) = rand(nq,1).*(max(X{k}) - min(X{k})) + min(X{k});
end
XQ = mat2cell(Xq,nq,ones(1,N));
[GRID{1:numel(X)}] = ndgrid(X{:});
clear out
for kk=1:nOut
    out(:,kk) = interpn(GRID{:},OUT{kk},XQ{:});
end
%%
obj = NPLR(N,X,OUT);
cplr = obj.eval(Xq);
RMSE(N) = sqrt(mean((out(:)-cplr(:)).^2))/sqrt(mean(out(:).^2));
disp(['RMSE = ' num2str(RMSE(N)) ' for N = ' num2str(N) '  in ' num2str(toc) ' (sec)'])
tic
semilogy(abs(out-cplr)+eps), hold on
end
grid on, box on, xlabel('q')
figure(2)
semilogy(RMSE,'.-b','MarkerSize',10),xlabel('Dimensions - N'), ylabel('RMSE'),grid on,box on