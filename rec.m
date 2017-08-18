function [g] = rec(hHat, Li, Lw, Lcit, tau, fs, t60)
% RMCLS with envelope constraint
%
%   g = rec(hHat, Li, Lw, Lcit, tau, fs, t60)
%
%   Input:
%       hHat      : estimated AIR
%       Li        : length of equalization filter
%       Lw        : length of relaxation window
%       Lcit      : sub-length of Lw to constrain to initial taps
%       tau       : delay
%       fs        : sampling frequency
%       t60       : reverberation time (secs)
%
%   Output:
%       g : equalizing filters
%
% Authors: F. Lim, W. Zhang
%
% Copyright (C) Imperial College London 2009-2010

[L, M] = size(hHat);

% Tolerance for error between the mask and the EIR
tol = 1e-3;    

% Convolution matrix
H = zeros(L+Li-1,M*Li);
for ii = 1:M
    H(:,(ii-1)*Li+1:ii*Li) = convmtx(hHat(:,ii),Li);
end
    
% Define mask according to Polack's method
if (Lcit > 0)
    n = Lcit:Lw-1;
else
    n = 1:Lw;
end
alpha = 10^-(3/t60/fs);
decurv_h = decacrv_v2(hHat(:,1));
beta = sqrt(decurv_h(Lw+1)*(1-alpha^2)/(alpha^(2*(tau+Lw))-alpha^(2*(tau+L))));
mask = beta*alpha.^n';
clear n;

% Define target impulse response  
if (Lcit > 0)
    hHatInitTaps = hHat(tau+1:tau+Lcit,1);
    d = [zeros(tau,1); hHatInitTaps; zeros(L+Li-Lcit-tau-1,1)];
else
    d = [zeros(tau,1); 1; zeros(L+Li-tau-2,1)];
end

% Define the weighting function
if (Lcit > 0)
    w = [ones(tau,1); ones(Lcit,1); zeros(Lw-Lcit,1); ones(L+Li-1-Lw-tau,1)];
else 
    % place one weight arbitrarily at the beginning to avoid the trivial solution.
    w = [ones(tau,1); 1; zeros(Lw-1,1); ones(L+Li-1-Lw-tau,1)];
end

if (Lcit > 0)
    nr_weights_before = sum(w);
else
    nr_weights_before = 0;
end
    
iterCount = 0;
while (true)
    iterCount = iterCount + 1;
    display(sprintf(['Iteration: ',num2str(iterCount), ...
                     ', Number of weights: ', num2str(sum(w))]));

    if (Lcit == 0 && iterCount == 2)
        w(tau+1) = 0;   % reset to zero the initial weight placed to avoid the trivial solution
    end
    g = rmcls(H,w,d);

    eir = zeros(length(hHat)+Li-1,1);
    for m = 1:M
        eir = eir + conv(hHat(:,m), g(:,m));
    end
    % Shape only the taps from Lcit+1:Lw
    eir = eir(tau+Lcit+1:tau+Lw);

    % update w
    wUpdated = false;
    while(sum(abs(eir)) > 0)
        [~,ind] = max(abs(eir));
        if (abs(eir(ind)) - mask(ind) > tol)
            % when updating w and d, compensate for the tau+Lcit taps not
            % considered in the comparison against the mask
            if (w(ind+tau+Lcit) == 0)
                w(ind+tau+Lcit) = 1;
                wUpdated = true;
                si = sign(eir(ind));
                d(ind+tau+Lcit) = si.*mask(ind);
            end
        end
        eir(ind) = 0;
    end
    
    
    if (~wUpdated)
        display(sprintf('w not updated'));
        break;
    end
end

nr_weights_after = sum(w);
total_nr_weights = nr_weights_after - nr_weights_before;
total_unconstr_weights = length(w) - sum(w);

function g = rmcls(H,w,d)
    W = repmat(w,1,M*Li);
    iH = pinv(W.*H);
    g = reshape(iH*(w.*d), Li, M); 
end

function decurv = decacrv_v2(h)
% decay curve computing
[Ld,Md] = size(h);
bsmatrx = triu(ones(Ld,Ld),0);
decurv = zeros(Ld, Md);
for md = 1:Md
    h(:,md) = h(:,md)/h(1,1);
    decurv(:,md) = bsmatrx * h(:,md).^2;
end
end

end