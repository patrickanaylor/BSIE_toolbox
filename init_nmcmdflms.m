function [shMDF H_hat] = init_nmcmdflms(L, Ms, Mch, xm0, H_hat0)

% This function initializes NMCMDFLMS algorithm
%
%	[shMDF H_hat] = init_nmcmdflms(L, Ms, M, H_hat0)
%
%	Input Parameters [size]:
%       L        : total filter length
%       Ms       : number of MDF segments
%       Mch      : number of channels
%       xm0      : input data
%       H_hat0   : initial filter coef. matrix (unit-norm constrained) [L x M]
%
%   Output parameters [size]:
%       shMDF    : MDF data structure
%       H_hat    : initialized filter coef. matrix [L x M]
%
% Authors: B. Castro and E.A.P. Habets
%
% History: 2011-03-01 Initial version
%
% Copyright (C) Bar-Ilan University 2011

error(nargchk(2,5,nargin));

if rem(2*L,Ms) ~= 0
   error('Segment length (2*L/Ms) needs to be an integer number.');
end
Ns = 2*L/Ms; % segment length

if nargin > 5 && ~isempty(H_hat0)
    if size(H_hat0,1) == L && size(H_hat0,2) == Mch
        H_hat = H_hat0;
    else
        error('hhat0 must be of size L times Mch.');
    end
else
    % Unit-norm initialization
    H_hat = [ones(1,Mch) ; zeros(L-1,Mch)]/sqrt(Mch);
end;

H1 = zeros(Ns , Mch , Ms);
for nn = 1:Mch
    hh1 = reshape(H_hat(:,nn),Ns/2,Ms); % rearange H_hat in blocks
    H1(:,nn,:) = fft(hh1,Ns);% padding each block
end

% Define shMDF data structure
shMDF.Ns = Ns;
shMDF.Ms = Ms;
shMDF.Mch = Mch;
shMDF.H_hat = H1;
shMDF.Xm = zeros(Ns, Mch, Ms); % Data buffer with overlap

% Initialization of shMDF.St
Xm = fft(xm0,shMDF.Ns);
P_x = conj(Xm).*Xm;
P_k = sum(P_x,2)*ones(1,Mch)-P_x;
shMDF.P_k_avg = P_k;