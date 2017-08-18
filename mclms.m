function [h_hat, J] = mclms(xin, h_hat, mu, ss_cntr)

% A time-domain Multi-channel LMS for blind channel identification [1]
%
%	[h_hat, J] = mclms(xin, h_hat, mu, ss_cntr)
%
%	Input Parameters [size]:
%       xin     : input vector [N x M]
%       h_hat   : current filter coef. matrix [L x M]
%       mu      : step-size
%       ss_cntr : type of step-size control (optional: default 'normalized'):
%                 'fixed' - fixed with unit-norm-constrained;
%                 'vss-unconstrained' - variable step-size;
%                 'vss' - variable step-size
%
%	Output parameters [size]:
%       h_hat   : updated filter coef. matrix [L x M]
%   	J       : value of cost function
%
%   References:
%       [1] Y. Huang and J. Benesty, "Adaptive multi-channel mean square and
%           Newton algorithms for blind channel identification", Signal Process.,
%           vol. 83, no. 8, pp. 1127-1138, Aug 2002
%
% Authors:  N.D. Gaubitch, E.A.P. Habets
%
% History:  2004-11-02 - Initial version by NG
%           2009-07-11 - Added various step-size control types by EH
%           2009-08-24 - For real and complex input vectors by EH
%                        Does not work with vss-normlized and vss!
%
% Copyright (C) Imperial College London 2004-2010

error(nargchk(3,4,nargin));

if nargin < 4
    ss_cntr = 'fixed';
end

e = xin.'*h_hat;  % cross-relation error 
e = e-e.';

% Note that norm(h_hat,2).^2 = 1
J = 0.5*(e(:)'*e(:));  % cost function
        
switch lower(ss_cntr)
    case 'fixed'
        DJ = 2*(conj(xin)*e-J*h_hat);
        h_hat = h_hat - mu*DJ;  % update
        h_hat = h_hat/norm(h_hat(:));  % constrain
        
    case 'vss-unconstrained'
        DJ = 2*(conj(xin)*e);
        mu_optimal = h_hat(:).' * DJ(:) / norm(DJ(:),2).^2;
        h_hat = h_hat - mu_optimal*DJ;  % update
        
    case 'vss'
        DJ = 2*(conj(xin)*e);
        mu_optimal = h_hat(:).' * DJ(:) / norm(DJ(:),2).^2;
        h_hat = h_hat - mu_optimal*DJ;  % update
        h_hat = h_hat/norm(h_hat(:));  % constrain
        
    otherwise
        error('Unknown step-size control type.');    
end       