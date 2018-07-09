function [g] = channel_shortening(h_hat, Li, Lw, k)

% Channel shortening for acoustic system equalization
%
%	[g] = channel_shortening(h_hat, Li, Lw, k)
%
%	Input Parameters [size]:
%       h_hat : M impulse responses of length L [L x M] 
%       Li    : length of the equalization filters
%       Lw    : length of the shortened impulse response
%       k     : delay of the target response
%
%   Output parameters [size]:
%       g     : equalization filters [Li x M] 
%
%   Reference:
%       [1] 
%
% Authors: W. Zhang
%
% History: 2009-07-06 - Initial version by W. Zhang
%
% Copyright (C) Imperial College London 2009-2010

% Initialization
[L M]= size(h_hat);
H = zeros(L+Li-1,M*Li); 
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(h_hat(:,ii),Li);
end

% Compute equalization system
wd = [zeros(k,1); ones(Lw,1); zeros(L+Li-Lw-k-1,1)];
wu = [ones(k,1); zeros(Lw,1); ones(L+Li-Lw-k-1,1)];
Wd = repmat(wd, 1, M*Li);
Wu = repmat(wu, 1, M*Li);
Bd = Wd.*H;
Au = Wu.*H;
clear wd wu Wd Wu;

B = Bd'*Bd;
A = Au'*Au;
clear Bd Au

[V,D] = eig(B, A);
ev = diag(D);
clear D;

% Find solution that produces minimum 2-norm
if(any(isinf(ev)))
    Vinf=V(:,isinf(ev));
    evall=H*Vinf;
    [~,I]=min(sqrt(sum(evall.^2,1)));
    gv=Vinf(:,I);
else
    [mev indmev] = max(ev);
    gv = V(:, indmev);
end
g = reshape(gv, Li, M);