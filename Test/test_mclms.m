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
% Authors: S. Lin, E.A.P. Habets
%
% History: 2009-06-25 Initial version by SL
%          2009-07-10 Major changes by EH
%
% Copyright (C) Imperial College London 2009-2010

clc;
clear; 
close all;

%% Initialization
M = 5;      % number of channels
fs = 8e3;   % sampling frequency
L = 128;    % channel length
SNR = 30;   % signal to noise ratio (dB)
N = 4*fs;   % data length (samples)
s_seed = 1; % seed for generating source signal
v_seed = 50;% seed for generating noise signal
mu = .8;    % step-size

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(M, L, fs, air, N, SNR, s_seed, v_seed);

% Initialize MCLMS
[xin, h_hat] = init_mclms(L, M);
npm_dB = zeros(N,1);
J = zeros(N,1);

%% Processing Loop: run MCLMS [2]
wbar = waitbar(0,'MCLMS');
ss_cntr = {'fixed','vss'};
h_hat = repmat(h_hat,[1 1 length(ss_cntr)]);
for nn = 1 : N
    waitbar(nn/N);
    xin = [x(nn, :); xin(1 : L-1, :)];
    for tt = 1 : length(ss_cntr)
        [h_hat(:,:,tt), J(nn,tt)] = mclms(xin, h_hat(:,:,tt), mu, ss_cntr{tt});
        npm_dB(nn,tt) = 20*log10(npm(h, h_hat(:,:,tt)));
    end
end
close(wbar);

%% Plot results
figure(1); 
phandle = plot_npm(npm_dB, fs, 1, '-o', 'MarkerIndices',1:floor(length(npm_dB)/20):length(npm_dB));
legend(['MCLMS' strcat(ss_cntr(2:end),'-MCLMS')]);
title(['L= ',num2str(L), ', \mu= ',num2str(mu), ...
    ', SNR= ',num2str(SNR), ', M= ', num2str(M)]);

figure(2); plot_J(J, fs);  % cost function
legend(['MCLMS' strcat(ss_cntr(2:end),'-MCLMS')]);
title(['L= ',num2str(L), ', \mu= ',num2str(mu), ...
    ', SNR= ',num2str(SNR), ', M= ', num2str(M)]);

figure(3); plot_filter(h, h_hat(:,:,2)) % filter coeff.