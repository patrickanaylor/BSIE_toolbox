% Run adaptive BSI algorithm: MCFLMS
% This script runs the MCFLMS algorithm for 
% blind system identification
%
% References:
%   [1] J. Allen and D. Berkley, "Image method for 
%       efficiently simulating small room acoustics",
%       JASA, 1979.
%   [2] Y. Huang and J. Benesty, "Frequency-Domain 
%       adaptive approaches to blind multi-channel
%       identification," IEEE Trans. Sig. Process., 
%       Vol. 51, No. 1, Jan 2003.
%
% Authors: E.A.P. Habets
%
% History: 2009-07-10 Initial version
%
% Copyright (C) Imperial College London 2009-2010

clc;
clear; 
close all;

%% Initialization
M = 5;     % number of channels
fs = 8e3;  % sampling frequency
L = 128;   % channel length
F = 2*L;   % frame lenght
SNR = 40;  % signal to noise ratio (in dB)
N = 5*fs;  % data length (samples)
s_seed = 1; % seed for generating source signal
v_seed = 50;% seed for generating noise signal
mu = 0.2;  % step-size

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(M, L, fs, air, N, SNR, s_seed, v_seed);

% Initiailize MCLMS
[h_hat] = init_mcflms(L, M);
ns = F-L+1;  % number of new samples per iteration
B = fix(N/ns);  % number of input blocks of length L
npm_dB = zeros(B,1);

%% Processing Loop: run MCFLMS [2]
wbar = waitbar(0,'MCFLMS');
ss_cntr = {'fixed','vss','rvss'};
h_hat = repmat(h_hat,[1 1 length(ss_cntr)]);
for bb = 1 : B
    waitbar(bb/B);
    if bb == 1
        xm = [zeros(F-ns,M); x(1:ns,:)]; 
    else
        xm = [xm(ns+1:end,:); x((bb-1)*ns+1:bb*ns,:)]; 
    end
    for tt = 1 : length(ss_cntr)
        [h_hat(:,:,tt)] = mcflms(xm, h_hat(:,:,tt), mu, ss_cntr{tt});
        npm_dB(bb,tt) = 20*log10(npm(h, h_hat(:,:,tt)));
    end
end
close(wbar);

%% Plot results
% NPM
figure(1); 
phandle = plot_npm(npm_dB, fs, ns);
lhandle = addmarkers(phandle,20);
legend(lhandle,['MCFLMS' strcat(ss_cntr(2:end),'-MCFLMS')]);
title(['L= ',num2str(L), ', \mu= ',num2str(mu), ...
    ', SNR= ',num2str(SNR), ', M= ', num2str(M)]);

% Time-domain coefficients
figure(2); plot_filter(h, h_hat(:,:,1));