% Run adaptive BSI algorithm: MCLMS
% This script runs the MCLMS algorithm for blind 
% system identification
%
% References:
%   [1] J. Allen and D. Berkley, "Image method for 
%       efficiently simulating small room acoustics",
%       JASA, 1979.
%   [2] Y. Huang and J. Benesty, "Adaptive multi-
%       channel mean square and Newton algorithms
%       for blind channel identification",
%       Signal Process., vol. 83, no. 8,
%       pp. 1127-1138, Aug 2002.
%
% Authors: E.A.P. Habets
%
% History: 2009-08-24 Initial version by EH
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: test_mclms_cmpl.m 425 2011-08-12 09:15:01Z mrt102 $

clc;
clear; 
close all;

%% Initialization
M = 3;      % number of channels
fs = 8e3;   % sampling frequency
L = 1;      % channel length
SNR = 40;   % signal to noise ratio (dB)
N = 0.2*fs; % data length (samples)
mu = .8;    % step-size

% Generate sensor signals and AIRs
s = randn(1,N) + 1i*randn(1,N);
h(:,1) = 1 * exp(1i*pi/2);
h(:,2) = 1 * exp(1i*0);
h(:,3) = 1 * exp(1i*2*pi/3);
for m = 1 : M
    x(:,m) = filter(h(:,m),1,s);
end

% Initiailize MCLMS
[xin, h_hat] = init_mclms(L, M);
npm_dB = zeros(N,1);
J = zeros(N,1);

%% Processing Loop: run MCLMS [2]
wbar = waitbar(0,'MCLMS');
%ss_cntr = {'fixed','normalized','vss'};
ss_cntr = {'normalized'};
h_hat = repmat(h_hat,[1 1 length(ss_cntr)]);
for nn = 1 : N
    waitbar(nn/N);
    xin = [x(nn, :); xin(1 : L-1, :)];
    for tt = 1 : length(ss_cntr)
        [h_hat(:,:,tt), J(nn,tt)] = mclms(xin, h_hat(:,:,tt), mu, ss_cntr{tt});
        %h_hat(:,:,tt) = real(h_hat_tmp);
        npm_dB(nn,tt) = 20*log10(npm(h, h_hat(:,:,tt)));
    end
end
close(wbar);

%% Plot results
figure(1); 
phandle = plot_npm(npm_dB, fs,1);
lhandle = addmarkers(phandle,20);
legend(lhandle,['MCLMS' strcat(ss_cntr(2:end),'-MCLMS')]);
title(['L= ',num2str(L), ', \mu= ',num2str(mu), ...
    ', SNR= ',num2str(SNR), ', M= ', num2str(M)]);

% figure(2); plot_J(J, fs);  % cost function
% legend(['MCLMS' strcat(ss_cntr(2:end),'-MCLMS')]);
% title(['L= ',num2str(L), ', \mu= ',num2str(mu), ...
%     ', SNR= ',num2str(SNR), ', M= ', num2str(M)]);

%figure(3); plot_filter(h, h_hat(:,:,1)) % filter coeff.

display([h ; h_hat]);
display(h_hat ./ h);