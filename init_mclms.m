function [xin, h_hat] = init_mclms(L, M, xin0, h_hat0)

% This function initializes MCLMS algorithm
%
%	[h_hat, xin] = init_mclms(L, N, M, h_hat0, xin0)
%
%	Input Parameters [size]:
%       L      : filter length
%       M      : number of channels
%       xin0   : initial input matrix [L x M]
%       h_hat0 : initial filter coef. matrix (unit-norm constrained) [L x M]
%
%   Output parameters [size]:
%       xin    : initialized signal matrix [L x M]
%       h_hat  : initialized filter coef. matrix [L x M]
%
% Authors: E.A.P. Habets
%
% History: 2009-07-10 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2009-2010

narginchk(2,4);

if nargin > 2 && ~isempty(xin0)
   if size(xin0,1) ~= L || size(xin0,2) ~= M
      error('xin0 must be of size N times M.');
   end
else			
   xin = zeros(L,M);
end;

if nargin > 3 && ~isempty(h_hat0)
   if size(h_hat0,1) ~= L || size(h_hat0,2) ~= M
      error('h_hat0 must be of size L times M.');
   end
else		
   % Unit-norm initialization
   h_hat = [ones(1,M) ; zeros(L-1,M)]/sqrt(M); 			
end;