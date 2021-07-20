clear; clc
N = 6;
S = @(p,q)sparse(NPLR.ShufflingMatrix(p,q));
Sv = @(p,q)(reshape(reshape(uint32(1:p*q),q,[])',[],1));
nSize = round(abs(rand(1,N))*2 + 3);
for k = 1:N
    n    = nSize(k);
    dx   = round(abs(rand*10 + 1));
    x0   = round(abs(rand*10 + 1));
    xf   = x0+dx*(n-1);
    X{k} = (x0:dx:xf)';
end
clear x0 xf dx n;
%%
temp_inc = X{1};
temp_dec = X{1};
for k=2:N
    temp_inc = kron(temp_inc,X{k});
    temp_dec = kron(X{k},temp_dec);
end
Z_inc = temp_inc;
Z_dec = temp_dec;
%% Complete S using either Sparse/Full Representation
tic
SS = S(prod(nSize(2:end)),nSize(1));
for k=2:N
temp = kron(eye(prod(nSize(1:k-1))) , S(prod(nSize(k+1:end)),nSize(k)));
SS = temp*SS;
end
toc
%% Only Vector of Indexes
tic
SSv = Sv(prod(nSize(2:end)),nSize(1));
for k=2:N-1
    m = prod(nSize(1:k-1));
    t_ =  Sv(prod(nSize(k+1:end)),nSize(k));
    pq = prod(nSize(k+1:end))*nSize(k);
%     tempv = t_;
%     for i=1:m-1
%         tempv = [tempv; t_ + i*pq];
%     end
    tempv = kron(uint32(ones(m,1)),t_) + uint32(kron((0:m-1)',pq*ones(pq,1)));
%     assert(all(tes_ == tempv))
    SSv = SSv(tempv);
end
toc
Ir = NPLR.IndexReversedOrder(nSize);
assert(all(SSv == Ir))
%%
norm(Z_inc - SS*Z_dec) 
norm(Z_inc - Z_dec(SSv)) 
whos SS*