% Run adaptive BSI algorithm: NMCFLMS_SC
% This script runs the NMCFLMS algorithm for
% blind system identification
%
% References:
%   [1] J. Allen and D. Berkley, "Image method for
%       efficiently simulating small room acoustics",
%       JASA, 1979.
%   [2] N. D. Gaubitch, Md. K. Hasan and P. A. Naylor, "Noise Robust
%       Adaptive Blind Channel Identification Using Spectral Constraints",
%       Proc. IEEE Int. Conf. on Acoust. Speech and Signal Processing
%       (ICASSP), Toulouse, France, May 2006.
%
% Authors: E.A.P. Habets
%
% History: 2010-03-27 Initial version
%
% Copyright (C) Imperial College London 2010

clc;
clear;
close all;

%% Initialization
M = 5;          % number of channels
fs = 8e3;       % sampling frequency
L = 128;        % channel length
F = 2*L;        % frame lenght
SNR = 20;       % signal to noise ratio (in dB)
N = 50*fs;      % data length (samples)
rho = 0.8;      % step-size
lambda = 0.8;   % forgetting-factor

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(M, L, fs, air, N, SNR);

x = [zeros(L,M); x];

% Initialize MCLMS_SC
[h_hat(:,:,1), P_k_avg(:,:,1)] = init_nmcflms(L, F, M, x(1:F,:));
[h_hat(:,:,2), P_k_avg(:,:,2), Pn] = init_nmcflms_sc(L, F, M, x(1:F,:));
ns = F-L+1;
B = fix(N/ns);  % number of input blocks of length ns
J = zeros(N,1);
npm_dB = zeros(B,1);

Xm = fft(x(1:F,:),F);
P_x = conj(Xm).*Xm;
delta = (M-1)*mean(mean(P_x));
%delta = sum(sum(P_x))/F*(M-1)/M/(10^(SNR/10));

%% Processing Loop: run NMCFLMS_SC
wbar = waitbar(0,'NMCFLMS with and without SC');
for bb = 1 : B
    waitbar(bb/B);
    xm = x((bb-1)*L+1:(bb+1)*L,:);
    
    [h_hat(:,:,1), P_k_avg(:,:,1), J_tmp] = nmcflms(xm, h_hat(:,:,1),...
        P_k_avg(:,:,1), rho, lambda, delta);
    J((bb-1)*ns+1:(bb)*ns,1) = J_tmp(F-L-ns+2:end);    
    npm_dB(bb,1) = 20*log10(npm(h, h_hat(:,:,1)));
    
     [h_hat(:,:,2), P_k_avg(:,:,2), Pn, J_tmp] = nmcflms_sc(xm, h_hat(:,:,2),...
         P_k_avg(:,:,2), Pn, rho, lambda, delta);
     J((bb-1)*ns+1:(bb)*ns,2) = J_tmp(F-L-ns+2:end);
     npm_dB(bb,2) = 20*log10(npm(h, h_hat(:,:,2)));
end
close(wbar);

%% Plot results
figure(1); plot_npm(npm_dB, fs, ns);  % NPM
legend('NMCFLMS','NMCFLMS-SC');
title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
    ', \lambda= ',num2str(lambda), ', SNR= ', ...
    num2str(SNR), ', M= ', num2str(M)]);

figure(2); plot_J(J, fs);  % cost function
legend('NMCFLMS','NMCFLMS-SC');
title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
    ', \lambda= ',num2str(lambda), ', SNR= ', ...
    num2str(SNR), ', M= ', num2str(M)]);

figure(3); plot_filter(h,h_hat(:,:,2));  % filter coeff.