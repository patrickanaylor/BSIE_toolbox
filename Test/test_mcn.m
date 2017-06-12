% Run adaptive BSI algorithm: MCN
% This script runs the MCN algorithm for blind
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
% History: 2009-07-10 Initial version
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: test_mcn.m 425 2011-08-12 09:15:01Z mrt102 $

clc;
clear; 
close all;

%% Initialization
M = 5;         % number of channels
fs = 8e3;      % sampling frequency
L = 16;        % channel length
SNR = 40;      % signal to noise ratio (in dB)
N = 0.2*fs;    % data length (samples)
s_seed = 1; % seed for generating source signal
v_seed = 50;% seed for generating noise signal
rho = 0.98;    % step-size
lambda = 0.99; % forgetting-factor

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(M, L, fs, air, N, SNR, s_seed, v_seed);

% Initiailize MCLMS
[xin, h_hat, R_hat] = init_mcn(L, M, [x(1,:) ; ...
    zeros(L-1,M)]);
npm_dB = zeros(N, 1);
N = size(x, 1);
J = zeros(N, 1);

%% Processing Loop: run MCN [2]
wbar = waitbar(0,'MCN');
for nn = 1 : N
    waitbar(nn/N);
    xin = [x(nn, :) ; xin(1 : L-1, :)];
    [h_hat, R_hat, J(nn)] = mcn(xin, h_hat,...
        R_hat, rho, lambda);
    npm_dB(nn) = 20*log10(npm(h, h_hat)); 
end
close(wbar);

%% Plot results
figure(1); plot_npm(npm_dB, fs);  % NPM
legend('MCN');
title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
    ', \lambda= ',num2str(lambda), ', SNR= ', ...
    num2str(SNR), ', M= ', num2str(M)]);

figure(2); plot_J(J, fs);  % cost function
title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
    ', \lambda= ',num2str(lambda), ', SNR= ', ...
    num2str(SNR), ', M= ', num2str(M)]);

figure(3); plot_filter(h,h_hat);  % filter coeff.
