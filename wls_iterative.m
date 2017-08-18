function [g, J] = wls_iterative(h_hat, Li, k, iter, w)

% Iterative weighted LS equalization system design using conjugate
% gradient method.
%
%	[g, J] = wls_iterative(h_hat, Li, k, iter, w)
%
%	Input Parameters [size]:
%       h_hat : M impulse responses of length L [L x M] 
%       Li    : length of the equalization filters
%       k     : delay of the target response
%       iter  : number of iterations
%       w     : weighting function [L+Li-1 x 1] (optional)
%
%   Output parameters [size]:
%       g     : equalization filters [Li x M] 
%       J     : cost function at each iteration
%
% Authors: W. Zhang
%
% History: 2009-07-06 - Initial version by W. Zhang
%
% Copyright (C) Imperial College London 2009-2010

% Initialization
[L, M] = size(h_hat);
H = zeros(L+Li-1,M*Li); 
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(h_hat(:,ii),Li);
end
d = [zeros(k,1); 1; zeros(L+Li-k-2,1)];

if nargin == 4
    w = ones(L+Li-1,1);
end
W = repmat(w,1,M*Li);
T = W.*H;
clear W;
F = T'*T;
md = T'*(w.*d);
clear T;

J = zeros(iter,1);

% Iterations
gv = zeros(M*Li,1); 
r0 = md - F*gv;
p1 = r0;
v = F*p1;
mu = r0'*r0/(p1'*v);
gv = gv + mu*p1;
r1 = r0 - mu*v;
J(1) = norm(w.*(H*gv-d))^2;
for j = 2:iter
    beta = r1'*r1/(r0'*r0);
    p2 = r1 + beta*p1;
    v = F*p2;
    mu = r1'*r1/(p2'*v);
    gv = gv + mu*p2;
    r0 = r1;
    r1 = r1 - mu*v;
    p1 = p2;
    J(j) = norm(w.*(H*gv-d))^2;
end
g = reshape(gv, Li, M);