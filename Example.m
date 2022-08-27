clearvars; close all;
load('ExampleData.mat');

X = {aoa,beta,cs4,cs5};
Y = {Cl,Cm,Cn};

%% Generate Test Data (within the limits)
% Random Data (Uniformally Distributed)
n = 1e4;
Xq = [rand(n,1)*(max(aoa)  - min(aoa))  + min(aoa),...
      rand(n,1)*(max(beta) - min(beta)) + min(beta),...
      rand(n,1)*(max(cs4)  - min(cs4))  + min(cs4),...
      rand(n,1)*(max(cs5)  - min(cs5))  + min(cs5)];
%% Evaluate Actual Data
XQ = mat2cell(Xq,length(Xq),ones(1,4));
[GRID{1:numel(X)}] = ndgrid(X{:});

for k=1:3
out(:,k) = interpn(GRID{:},Y{k},XQ{:});
end
%% CPLR
obj = PMLR(4,X,Y); % Create Object
Ypmlr = obj.eval(Xq);   % Evaluate at Query points

% Compute RMSE
RMSE_pmlr = rms(out-Ypmlr);  rRMSE_pmlr = RMSE_pmlr./rms(out);

%% Display
fprintf('PMLR\n');
fprintf(' RMSE_PMLR = [%2.6e   %2.6e   %2.6e]\n', RMSE_pmlr)
fprintf('rRMSE_PMLR = [%2.3e(%%)   %2.3e(%%)   %2.3e(%%)]\n',rRMSE_pmlr*100)
