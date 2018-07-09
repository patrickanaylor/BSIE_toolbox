function [g] = wls(h_hat, Li, k, w)

% Weighted LS equalization system design for both single and 
% multichannel acoustic systems. 
%
%	[g] = wls(h_hat, Li, k, w)
%
%	Input Parameters [size]:
%       h_hat : M impulse responses of length L [L x M] 
%       Li    : length of the equalization filters
%       k     : delay of the target response
%       w     : weighting function [L+Li-1 x 1] (optional)
%
%   Output parameters [size]:
%       g     : equalization filters [Li x M] 
%
%   Reference:
%       [1] M. Miyoshi .etc, "Inverse filtering of room acoustics", 
%           IEEE Trans. ASSP, vol. 36, 1988.
%
%       [2] Y.Huang, J. Benesty and J. Chen, "A Blind Channel Identification-Based
%           Two-Stage Approach to Separation and Dereverberation of Speech
%           Signals in a Reverberant Environment," IEEE Trans. Speech Audio
%           Processing, vol. 13, no. 5 pp. 882-895, 2005.
%
% Authors: W. Zhang
%
% History: 2009-07-06 - Initial version by W. Zhang
%
% Copyright (C) Imperial College London 2009-2010

[L, M] = size(h_hat);

% Define target impulse response
d = [zeros(k,1); 1; zeros(L+Li-k-2,1)];

% H is the matrix in Equation (11a,11b) in [1]
H = zeros(L+Li-1,M*Li);
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(h_hat(:,ii),Li);
end

% Compute inverse
if nargin == 3
    iH = pinv(H);
    g = reshape(iH*d, Li, M);
else
    W = repmat(w,1,M*Li);
    iH = pinv(W.*H);
    g = reshape(iH*(w.*d), Li, M);
end