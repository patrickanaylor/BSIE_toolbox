function [h_hat, P_k_avg, J] = nmcflms(xm, h_hat, P_k_avg, rho, lambda, delta)

% The blind Normalized Multichannel Frequency-domain LMS (NMCFLMS)
%
%   [h_hat, P_k_avg, J] = nmcflms(xm, h_hat, P_k_avg, rho, lambda, SNR)
%
%   Input Parameters [size]:
%       xm      : input matrix [F x M]
%       h_hat   : current filter coef. matrix [L x M]
%       P_k_avg : estimated PSD matrix [F x M]
%       rho     : step size (optional: default rho=0.8)
%       lambda  : forgetting factor (0 < lambda < 1) (optional)
%       delta   : regularization (optional)
%
%   Output Parameters:
%       h_hat   : updated filter coef. matrix [L x M]
%       P_k_avg : updated PSD matrix [F x M]
%       J       : cost function values [F-L+1 x 1]
%
%   Remarks:
%       This algorithm is suitable for estimating IRs generated using the
%       method of images.
%
%       For simulated IRs, leading zeros should be removed before employing
%       this algorithm, i.e., identify the first direct path and remove all
%       values before it and apply same process to other channels. In order
%       to keep IRs the same length, random coefficients with small values
%       can be appended at the tail part.
%
%    References:
%       [1] Y. Huang and J. Benesty, "Frequency-Domain adaptive approaches to
%           blind multi-channel identification," IEEE Trans. Sig. Process. Vol. 51
%           No. 1, Jan 2003.
%
%   Authors: N. D. Gaubitch, E.A.P. Habets
%
%   History: 2006-03-29 Initial version by N. D. Gaubitch
%            2009-07-11 Minor Changes by E.A.P. Habets
%
% Copyright (C) Imperial College London 2006-2010

narginchk(3,6);

% Initialization
L = size(h_hat,1);
F = size(xm,1);
M = size(h_hat,2);
U = zeros(F,M);
J = zeros(F-L+1,1);
e = zeros(F,M);  % cross-relation error

if nargin < 4 || isempty(rho)
    rho = 0.8;
end
if nargin < 5 || isempty(lambda)
    lambda  = (1 - (1/(3*L)))^L;
end

h_hat_10 = fft(h_hat,F);
Xm = fft(xm);
P_x = conj(Xm).*Xm;
P_k = sum(P_x,2)*ones(1,M)-P_x;
P_k_avg = lambda*P_k_avg + (1-lambda)*P_k;

if nargin < 6 || isempty(delta)
    delta = (M-1)*mean(mean(P_x));
    % delta = sum(sum(P_x))/F*(M-1)/M/(10^(SNR/10));  % regularization
end

% Calculate cross-relation error
for kk = 1:M
    for ll = 1:M
        e(:,ll) = Xm(:,ll).*h_hat_10(:,kk)-Xm(:,kk).*h_hat_10(:,ll);
    end
    e = real(ifft(e));
    e(1:L-1,:) = 0;
    J = J + 0.5*sum(abs(e(L:end,:)).^2, 2);
    e = fft(e);
    U(:,kk) = sum(conj(Xm).*e,2);
end

% Update h_hat
h_hat_10 = h_hat_10 - rho./(P_k_avg + delta).*U;
h_hat = real(ifft(h_hat_10));
h_hat = h_hat(1:L,:);
h_hat = h_hat/norm(h_hat(:));