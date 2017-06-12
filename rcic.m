function g = rcic(hHat, Li, Lw, Lws, tau, desiredEq)
% rcic RMCLS with constrained initial coefficients
% 
%   USAGE: 
%       g = rcic(hHat, Li, Lw, Lws, tau, desiredEq)
% 
%   INPUT:
%       hHat      : M impulse responses of length L [L x M] 
%       Li        : length of filter
%       Lw        : length to channel shorten the EIR to
%       Lws       : length of initial taps to constrain
%       tau       : delay of the target response
%       desiredEq : desired equalized response (optional)
%
%   OUTPUT:
%       g : equalizing filter
%
%   REFERENCES:
%
%   AUTHOR   :  Felicia Lim
%**************************************************************************

%% Initialization
if (nargin == 5)
    % Choose the first estimated channel by default
    desiredEq = hHat(:,1);
end

[L, M] = size(hHat);

%% Define the desired target response
d = [zeros(tau,1); desiredEq(1:Lws); zeros(L+Li-Lws-tau-1,1)];


%% Calculate the convolution matrix
H = zeros(L+Li-1,M*Li);
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(hHat(:,ii),Li);
end


%% Define the relaxation window
if (Lw > 0)
    w = [ones(tau,1); ones(Lws,1); zeros(Lw-Lws,1); ones(L+Li-1-Lw-tau,1)];
else
    w = ones(tau+L+Li-1,1);
end

W = repmat(w,1,M*Li);

%% Find the equalizing filters
iH = pinv(W.*H);
g = reshape(iH*(w.*d),Li,M);

end

