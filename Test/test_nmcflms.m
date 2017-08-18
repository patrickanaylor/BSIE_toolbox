% Run adaptive BSI algorithm: NMCFLMS
% This script runs the NMCFLMS algorithm for 
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
%close all;

%% Initialization
M = 5;          % number of channels
fs = 8e3;       % sampling frequency
L = 128;        % channel length
F = 2*L;        % frame lenght
SNR = 40;       % signal to noise ratio (in dB)
s_seed = 1;     % seed for generating source signal
v_seed = 50;    % seed for generating noise signal
N = 5*fs;       % data length (samples)
rho = 0.8;      % step-size
lambda = 0.98;  % forgetting-factor

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(M, L, fs, air, N, SNR, s_seed, v_seed);

% Initiailize NMCFLMS
[h_hat, P_k_avg] = init_nmcflms(L, F, M, x(1:F,:));
ns = F-L+1;
B = fix(N/ns);  % number of input blocks of length L
J = zeros(N,1);
npm_dB = zeros(B,1);

Xm = fft(x(1:F,:),F);
P_x = conj(Xm).*Xm;
delta = (M-1)*mean(mean(P_x));

%% Processing Loop: run NMCFLMS [2]
wbar = waitbar(0,'NMCFLMS');
for bb = 1 : B
    waitbar(bb/B);
    if bb == 1
        xm = [zeros(F-ns,M); x(1:ns,:)]; 
    else
        xm = [xm(ns+1:end,:); x((bb-1)*ns+1:bb*ns,:)]; 
    end
    [h_hat, P_k_avg, J_tmp] = nmcflms(xm, h_hat,...
        P_k_avg, rho, lambda, delta);
    J((bb-1)*ns+1:(bb)*ns) = J_tmp(F-L-ns+2:end);
    npm_dB(bb) = 20*log10(npm(h, h_hat));   
end
close(wbar);

%% Plot results
figure(1); plot_npm(npm_dB, fs, ns);  % NPM
legend('NMCFLMS');
title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
    ', \lambda= ',num2str(lambda), ', SNR= ', ...
    num2str(SNR), ', M= ', num2str(M)]);

figure(2); plot_J(J, fs);  % cost function
legend('NMCFLMS');
title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
    ', \lambda= ',num2str(lambda), ', SNR= ', ...
    num2str(SNR), ', M= ', num2str(M)]);

figure(3); plot_filter(h,h_hat);  % filter coeff.