function [h_hat, P_k_avg, Pn, J] = nmcflms_sc(xm, h_hat, P_k_avg, Pn, rho, lambda, delta)

% The blind Normalized Multichannel Frequency-domain LMS with Spectral Constraint (NMCFLMS_SC)
%
%   [h_hat, P_k_avg, J] = nmcflms(xm, h_hat, P_k_avg, rho, lambda, SNR)
%
%   Input Parameters [size]:
%       xm      : input matrix [F x M]
%       h_hat   : current filter coef. matrix [L x M]
%       P_k_avg : estimated PSD matrix [F x M]
%       Pn      : penality matrix [F x M]
%       rho     : step size (optional: default rho=0.8)
%       lambda  : exponential forgetting factor (0 < lambda < 1) (optional)
%       delta   : regularization (optional)
%
%   Output Parameters:
%       h_hat   : updated filter coef. matrix [L x M]
%       P_k_avg : updated PSD matrix [F x M]
%       Pn      : updated penality matrix [F x M]
%       J       : cost function values [F-L+1 x 1]
%
%   Remarks:
%       This algorithm is suitable for estimating IRs generated using the
%       method of images.
%
%       For simulated IRs, leading zeros should be removed before employing
%       this algorithm, i.e., identify the first direct path and remove all
%       values before it and apply same process to other channels.
%
%    References:
%       [1] N. D. Gaubitch, Md. K. Hasan and P. A. Naylor, "Noise Robust
%           Adaptive Blind Channel Identification Using Spectral
%           Constraints", Proc. IEEE Int. Conf. on Acoust. Speech and 
%           Signal Processing (ICASSP), Toulouse, France, May 2006.
%
%   Authors: N.D. Gaubitch and E.A.P. Habets
%
%   History: 2010-03-27 Initial version by N.D. Gaubitch and E.A.P. Habets
%
% Copyright (C) Imperial College London 2010
% Version: $Id: nmcflms_sc.m 425 2011-08-12 09:15:01Z mrt102 $

error(nargchk(4,7,nargin));

% Initialization
L = size(h_hat,1);
F = size(xm,1);
M = size(h_hat,2);
U = zeros(F,M);
J = zeros(F-L+1,1);
e = zeros(F,M);  % cross-relation error

if nargin < 5 || isempty(rho)
    rho = 0.8;
end
if nargin < 6 || isempty(lambda)
    lambda  = (1 - (1/(3*L)))^L;
end

h_hat_10 = fft(h_hat,F);
Xm = fft(xm);
P_x = conj(Xm).*Xm;
P_k = sum(P_x,2)*ones(1,M)-P_x;
P_k_avg = lambda*P_k_avg + (1-lambda)*P_k;

if nargin < 7 || isempty(delta)
    delta = (M-1)*mean(mean(P_x));
    % delta = sum(sum(P_x))/F*(M-1)/M/(10^(SNR/10));  % regularization
end

% Initialization spectral constraint
wt      = impz([1 -1],1,M*L);
WT      = abs(fft(wt)).^4;
WT      = WT/sum(WT);
gamma   = 0.5;
beta    = 0.0075;

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
h_hat_10 = h_hat_10 - rho * (1./(P_k_avg + delta).*U + Pn);
h_hat = real(ifft(h_hat_10));
h_hat = h_hat(1:L,:);
h_hat = h_hat/norm(h_hat(:));

% Update Pn
hhatv = h_hat(:); 
Hhatv = fft(hhatv,M*L);
E = hhatv.'*hhatv; % Note: Result is approximatly one due to unit-norm constraint
WH = WT.*Hhatv;
Gh = M*L*real(ifft(WH));
ep = 1/(M*L)*hhatv.'*Gh - gamma*E;  
dep = reshape(Gh,L,M);
Pn = 4*beta*conj(fft(ep*dep,F)); 