% Run subspace BSI algorithm: MCLS
% This script runs the MCLS algorithm for blind 
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
% History: 2010-03-23 Initial version by EH
%
% Copyright (C) Imperial College London 2010

clc;
clear; 
close all;

%% Initialization
M = 3;       % number of channels
fs = 8e3;    % sampling frequency
L = 3;       % channel length
N = 2*fs;    % data length (samples)
N_ittr = 10; % number of itterations
SNR = 30;

%% Generate impulse responses
% theta = pi/10; % Ill Conditioned System
% phi = pi/10; % Ill Conditioned System
theta = pi/10; % Well Conditioned System
phi = 2*pi/3; % Well Conditioned System
h(:,1) = [1 -2*cos(theta) 1]';
h(:,2) = [1 -2*cos(theta+phi) 1]';
h(:,3) = [1 -2*cos(theta+2*phi) 1]';

%% Generate signals
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));
s = randi(2,1,N);   % Generate N length of either 0 or 1.
s = s - mean(s);
for i = 1:M
    x(:,i) = filter(h(:,i),1,s);
end

%% Add noise
SNR_M = SNR + [0 0 -15];
for i = 1:M
    x(:,i) = awgn(x(:,i),SNR_M(i));
end

%% Processing Loop: run MCLS
[h_hat(:,:,1), R_hat(:,:,1)] = mcls(x,L,'recursive_alt');
npm_dB(1) = 20*log10(npm(h, h_hat(:,:,1)));

for ittr = 2 : N_ittr
    [h_hat(:,:,ittr),R_hat(:,:,ittr)] = mcls(x,L,'spcc',h_hat(:,:,ittr-1));
    npm_dB(ittr) = 20*log10(npm(h, h_hat(:,:,ittr)));
end

%% Plot results
figure(1);
plot(npm_dB,'b-o');
grid on;
xlabel('Itterations');
ylabel('NPM (dB)')