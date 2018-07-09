function [h_hat, P_k, Pn] = init_nmcflms_sc(L, F, M, xm0, h_hat0)

% This function initializes NMCFLMS_SC algorithm
%
%	[h_hat, P_k, Pn] = init_nmcflms_sc(L, F, M, xm0, h_hat0)
%
%	Input Parameters [size]:
%       L      : filter length
%       F      : frame length
%       M      : number of channels
%       xm0    : initial input block used to calculate P_k [F x M]
%       h_hat0 : initial filter coef. matrix (unit-norm constrained) [L x M]
%
%   Output parameters [size]:
%       h_hat  : initialized filter coef. matrix [L x M]
%       P_k    : initialized PSDs [F x M]
%       Pn     : initialized P_n [F x M]
%
% Authors: E.A.P. Habets
%
% History: 2010-03-27 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2009-2010

narginchk(3,5);

if nargin > 3 && ~isempty(xm0)
    if size(xm0,1) == F && size(xm0,2) == M
        Xm = fft(xm0,F);
        P_x = conj(Xm).*Xm;
        P_k = sum(P_x,2)*ones(1,M)-P_x;        
    else
        error('xin0 must be of size F times M.');
    end
else
    P_k = zeros(F,M);
end;

if nargin > 4 && ~isempty(h_hat0)
    if size(h_hat0,1) == L && size(h_hat0,2) == M
        h_hat = h_hat0;
    else
        error('h_hat0 must be of size L times M.');
    end
else
    % Unit-norm initialization
    h_hat = [ones(1,M) ; zeros(L-1,M)]/sqrt(M);
end;

Pn = zeros(F,M);