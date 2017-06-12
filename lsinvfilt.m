function [g] = lsinvfilt(h_hat, Li, k)

% LS equalization system design for both single and 
% multichannel acoustic systems. 
%
%	[g] = lsinvfilt(h_hat, Li, k)
%
%	Input Parameters [size]:
%       h_hat : M impulse responses of length L [L x M] 
%       Li    : length of the equalization filters
%       k     : delay of the target response
%
%   Output parameters [size]:
%       g     : equalization filters [Li x M] 
%
%	Remarks: 
%       When Li = L-1 we obtain the MINT solution as proposed in [1].
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
% Version: $Id: lsinvfilt.m 425 2011-08-12 09:15:01Z mrt102 $

[L, M] = size(h_hat);

% Define target impulse response
d = [zeros(k,1); 1; zeros(L+Li-k-2,1)];

% H is the matrix in Equation (11a,11b) in [1]
%fprintf('Generating convmtx\n');
%tic;
H = zeros(L+Li-1,M*Li);
%H2=[];
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(h_hat(:,ii),Li);
    %H2=[H2 myconvmtx(h_hat(:,ii),Li)];
end
%toc;

% Compute inverse

% if(size(H,1)==size(H,2))
% %     fprintf('Computing MINT inverse LU');
% %     tic;
% %     [F.L F.U p]=lu(H,'vector');     % Slow bit (need to exploit block toeplitz)
% %     F.p=sparse(1:length(p),p,1);
% %     opL.LT=true;
% %     opU.UT=true;
% %     g=reshape(linsolve(F.U,linsolve(F.L,F.p*d,opL),opU),Li,M);
% %     g=reshape(inverse(H)*d,Li,M);
% %     toc
%     
%     fprintf('Computing MINT inverse with \\ \n');
%     %tic
%     g = reshape(H\d, Li, M);
%     %toc
%     
%     %max(max(g2-g))
% else
%     fprintf('Computing MINT inverse with pinv(H) \n');
%     
%     tic
    iH=pinv(H);
    g = reshape(iH*d, Li, M);
%     toc
% end
%
% function M = myconvmtx(x,nh)
% n = length(x);
% M = sparse(bsxfun(@plus,(1:n)',0:(nh-1)), ...
%     repmat(1:nh,n,1),repmat(x,1,nh),n+nh-1,nh);
