function [h_hat, J] = mcflms(xm, h_hat, mu, ss_cntr)

% The blind Multichannel Frequency-domain LMS (MCFLMS)
%
%   [h_hat, J] = mcflms(xm, h_hat, mu, ss_cntr)
%
%   Input Parameters [size]:
%       xm      : input matrix [F x M]
%       h_hat   : current filter coef. matrix [L x M]
%       mu      : step size (optional: default mu=0.9)
%       ss_cntr : type of step-size control (optional: default 'fixed'):
%                 'fixed' - fixed with unit-norm-constrained;
%                 'vss' - variable step-size
%                 'rvss' - robust variable step-size
%
%   Output Parameters:
%       h_hat   : updated filter coef. matrix [L x M]
%       J       : value of the cost function
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

narginchk(2,4);

% Initialization
L = size(h_hat,1);
F = size(xm,1);
M = size(h_hat,2);
U = zeros(F,M);
J = zeros(F-L+1,1);
e = zeros(F,M);  % cross-relation error

if nargin < 3 || isempty(mu)
    mu = 0.9;
end
if nargin < 4 
    ss_cntr = 'fixed';
end

h_hat_10 = fft(h_hat,F);
Xm = fft(xm);

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

% Update and constrain h_hat
switch lower(ss_cntr)
    case 'fixed'
        h_hat_1x = h_hat_10 - mu*U;
        h_hat = real(ifft(h_hat_1x));
        h_hat = h_hat(1:L,:);
        h_hat = h_hat/norm(h_hat(:));
    
    case 'vss'
        mu_optimal = h_hat_10(:)' * U(:) ./ norm(U(:),2).^2;
        h_hat_1x = h_hat_10 - mu_optimal*U;
        h_hat = real(ifft(h_hat_1x));
        h_hat = h_hat(1:L,:);
        h_hat = h_hat/norm(h_hat(:));
    
    case 'rvss'
        DJ_p = 2./abs(h_hat_10).^2 .* h_hat_10;
        beta = abs(DJ_p(:)' * U(:) ./ norm(DJ_p(:),2).^2);
        mu_optimal = h_hat_10(:)' * U(:) ./ norm(U(:),2).^2;
        h_hat_1x = h_hat_10 - mu_optimal*U + mu_optimal*beta*DJ_p;  % update
        h_hat = real(ifft(h_hat_1x));
        h_hat = h_hat(1:L,:);    
        h_hat = h_hat/norm(h_hat(:));
        
    otherwise
        error('Unknown step-size control type.');
end