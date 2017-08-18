function [g] = rmcls(h_hat, Li, Lw, k)

% RMCLS equalization system design. 
%
%	[g] = wls(h_hat, Li, Lw, k)
%
%	Input Parameters [size]:
%       h_hat : M impulse responses of length L [L x M] 
%       Li    : length of the equalization filters
%       Lw    : length of relaxed window
%       k     : delay of the target response
%
%   Output parameters [size]:
%       g     : equalization filters [Li x M] 
%
% Authors: W. Zhang
%
% Copyright (C) Imperial College London 2009-2010

[L, M] = size(h_hat);

% Define target impulse response
d = [zeros(k,1); 1; zeros(L+Li-k-2,1)];

w = [ones(k,1); 1; zeros(Lw-1,1); ones(L+Li-1-Lw-k,1)];

% H is the matrix in Equation (11a,11b) in [1]
H = zeros(L+Li-1,M*Li);
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(h_hat(:,ii),Li);
end

W = repmat(w,1,M*Li);

iH = pinv(W.*H);
g = reshape(iH*(w.*d), Li, M);

% tic
% % Inspired by Factorize
% [F.L F.U p]=lu(W.*H,'vector');     % Slow bit (need to exploit block toeplitz)
% F.p=sparse(1:length(p),p,1);
% opL.LT=true;
% opU.UT=true;
% g=reshape(linsolve(F.U,linsolve(F.L,F.p*(w.*d),opL),opU),Li,M);
% toc