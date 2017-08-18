function [h_hat, P_k_avg, J] = rnmcflms(xm, h_hat, P_k_avg, rho, lambda, delta)

% The blind Robust Normalized Multichannel Frequency-domain LMS (RNMCFLMS)
%
%   [h_hat, P_k_avg, J] = nmcflms(xm, h_hat, P_k_avg, rho, lambda, SNR)
%
%   Input Parameters [size]:
%       xm      : input matrix [F x M]
%       h_hat   : current filter coef. matrix [L x M]
%       P_k_avg : estimated PSD matrix [F x M]
%       rho     : step size (optional: default rho=0.8)
%       lambda  : exponential forgetting factor 
%                 (0 < lambda < 1) (optional)
%       delta   : regularization (optional)
%
%   Output Parameters:
%       h_hat   : updated filter coef. matrix [L x M]
%       P_k_avg : updated PSD matrix [F x M]
%       J       : cost function values [F-L+1 x 1]
%
%    References:
%       [1] Y. Huang and J. Benesty, "Frequency-Domain adaptive approaches to
%           blind multi-channel identification," IEEE Trans. Sig. Process. 
%           vol. 51 no. 1, Jan 2003.
%       [2] M. Haque and M. Hasan, "Noise robust multichannel 
%           frequency-domain LMS algorithms for blind channel 
%           identification," IEEE Signal Process. Lett., vol. 15,
%           pp. 305-308, 2008.
%
%   Authors: E.A.P. Habets
%
%   History: 2009-08-05 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2009-2010

error(nargchk(3,6,nargin));

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
    delta = (M-1) * mean(mean(P_x));
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
DJ_p = 2./abs(h_hat_10).^2 .* h_hat_10;
DJ_f = U./(P_k_avg+delta);
beta = abs(DJ_p(:)' * DJ_f(:) ./ norm(DJ_p(:),2).^2);
h_hat_10 = h_hat_10 - rho*DJ_f + rho*beta*DJ_p;
h_hat = real(ifft(h_hat_10));
h_hat = h_hat(1:L,:);
h_hat = h_hat/norm(h_hat(:));