function [h_hat] = init_mcflms(L, M, h_hat0)

% This function initializes MCFLMS algorithm
%
%	[h_hat] = init_mcflms(L, M, h_hat0)
%
%	Input Parameters [size]:
%       L      : filter length
%       M      : number of channels
%       h_hat0 : initial filter coef. vector (unit-norm constrained) [L x M]
%
%   Output parameters [size]:
%       h_hat  : initialized filter coef. matrix [L x M]
%
% Authors: E.A.P. Habets
%
% History: 2009-07-10 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2009-2010

narginchk(2,3);

if nargin > 2 && ~isempty(h_hat0)
    if size(h_hat0,1) == L && size(h_hat0,2) == M
        h_hat = h_hat0;
    else
        error('hhat0 must be of size L times M.');
    end
else
    % Unit-norm initialization
    h_hat = [ones(1,M) ; zeros(L-1,M)]/sqrt(M);
end;