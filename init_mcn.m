function [xin, h_hat, Rhat] = init_mcn(L, M, xin0, h_hat0)

% This function initializes MCLMS algorithm
%
%	[xin, h_hat, Rhat] = init_mcn(L, N, M, xin0, h_hat0)
%
%	Input Parameters [size]:
%       L      : filter length
%       M      : number of channels
%       xin0   : initial input matrix [L x M]
%       h_hat0 : initial filter coef. matrix (unit-norm constrained) [L x M]
%
%   Output parameters [size]:
%       xin    : initialized input matrix [L x M]
%       h_hat  : initialized filter coef. matrix [L x M]
%       Rhat   : covariance matrix [M L x M L]
%
% Authors: E.A.P. Habets
%
% History: 2009-07-10 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: init_mcn.m 425 2011-08-12 09:15:01Z mrt102 $

error(nargchk(2,4,nargin));

if nargin > 2 && ~isempty(xin0)
    if size(xin0,1) == L && size(xin0,2) == M
        xin = [xin0(1,:); zeros(L-1,M)];
        Rhat = trace(xin0'*xin0)/L * eye(M*L);
    else
        error('xin0 must be of size L times M.');
    end
else
    xin = zeros(L,M);
    Rhat = eye(M*L);
end;

if nargin > 3 && ~isempty(h_hat0)
    if size(h_hat0,1) ~= L || size(h_hat0,2) ~= M
        error('h_hat0 must be of size L times M.');
    end
    h_hat = h_hat0;
else
    % Unit-norm initialization
    h_hat = [ones(1,M) ; zeros(L-1,M)]/sqrt(M);
end;