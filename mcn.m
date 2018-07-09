function [h_hat, R_hat, J] = mcn(xin, h_hat, R_hat, rho, lamda)

% Multi-channel Newton algorithm for blind channel identification [1]
%
%	[h_hat, R_hat, J] = mcn(xin, h_hat, R_hat, rho, lamda)
%
%   Input parameters:
%   	xin   : input matrix [L x M]
%       h_hat : current filter coef. matrix [L x M]
%       R_hat : covariance matrix [M L x M L]
%       rho   : step size (optional: default rho=0.95)
%       lamda : exponential forgetting factor (0 < lamda < 1)
%               (optional: default lamda=0.99)
%
%   Outputs parameters:
%       h_hat : updated filter coef. matrix [L x M]
%       R_hat : covariance matrix [M L x M L]
%       J     : value of cost function
%
%   References:
%       [1] Y. Huang and J. Benesty, "Adaptive multi-channel mean square and
%           Newton algorithms for blind channel identification", Signal Process.,
%           vol. 83, no. 8, pp. 1127-1138, Aug 2002
%
% Authors: N.D. Gaubitch, E.A.P. Habets
%
% History: 2003-05-15 Initial Version by N.D. Gaubitch
%        : 2004-10-18 Minor changes by N.D. Gaubitch
%        : 2009-07-10 Minor changes by E.A.P. Habets
%
% Copyright (C) Imperial College London 2003-2010

narginchk(3,5);

L = size(h_hat,1);
M = size(xin,2);
h_hat_ud = h_hat(:);
if nargin < 4   
    rho = 0.95;
end                     
if nargin < 5
    lamda = 0.99;
end                     
epsilon = 10^-8;

R_tilde = -xin(:)*xin(:)';
Rr = xin*xin';
for jj = 1:M
   row = (jj-1)*L+1:jj*L;
   for ii = 1:M
       col = (ii-1)*L+1:ii*L;
       R_tilde(row,col) = R_tilde(row,col)';
   end
   R_tilde(row,row) = R_tilde(row,row) + Rr;
end        
R_hat = lamda*R_hat + R_tilde;

H_mat = h_hat_ud*h_hat_ud';
W = R_hat - 2*H_mat*R_hat - 2*R_hat*H_mat;
B = inv(W+epsilon*eye(M*L,M*L)./(1+epsilon));

E = xin'*h_hat;
e = triu(E-E'); 
chi = e(:)'*e(:);
h_hat_ud = h_hat_ud - rho*B*(R_tilde*h_hat_ud - chi*h_hat_ud);
h_hat_ud = h_hat_ud/norm(h_hat_ud);   
J = chi/(h_hat_ud'*h_hat_ud);
h_hat = reshape(h_hat_ud,L,M);