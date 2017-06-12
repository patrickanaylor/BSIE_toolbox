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
% History: 2009-09-02 Initial version by EH
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: test_mclms_simple.m 425 2011-08-12 09:15:01Z mrt102 $

clc;
clear; 
close all;

%% Initialization
M = 3;              % number of channels
fs = 8e3;           % sampling frequency
L = 3;              % channel length
SNR = [30 30 20];   % signal to noise ratio (dB)
N = .5*fs;          % data length (samples)
mu = .4;            % step-size

%% Generate simple impulse responses
theta = pi/10;
phi = pi/10; % IC
%phi = 2*pi/3; % WC

h(:,1) = [1 -2*cos(theta) 1]';
h(:,2) = [1 -2*cos(theta+phi) 1]';
h(:,3) = [1 -2*cos(theta+2*phi) 1]';

%% Generate signals
s = randint(1,N);   % Generate D length of either 0 or 1.
s = s - mean(s);
for i = 1:M
    x(:,i) = filter(h(:,i),1,s);
    x(:,i) = awgn(x(:,i),SNR(i));
end

%% Initiailize MCLMS
global R_hat;
R_hat = zeros(L*M,L*M);
for tf=1:N-L+1
    xin = x(tf:tf+L-1,:);
    R_hat = 0.98 * R_hat + (1-0.98) * xin(:)*xin(:)';
end
[xin, h_hat] = init_mclms(L, M);
npm_dB = zeros(N,1);
J = zeros(N,1);

%% Processing Loop: run MCLMS [2]
wbar = waitbar(0,'MCLMS');
ss_cntr = {'fixed','normalized','vss'};
h_hat = repmat(h_hat,[1 1 length(ss_cntr)]);
for nn = 1 : N
    waitbar(nn/N);
    xin = [x(nn, :); xin(1 : L-1, :)];
    for tt = 1 : length(ss_cntr)
        if tt == 1
            [h_hat(:,:,tt), J(nn,tt)] = mclms(xin, h_hat(:,:,tt), mu/10, ss_cntr{tt});
        else
            [h_hat(:,:,tt), J(nn,tt)] = mclms(xin, h_hat(:,:,tt), mu, ss_cntr{tt});
        end
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

figure(2); plot_J(J, fs);  % cost function
legend(['MCLMS' strcat(ss_cntr(2:end),'-MCLMS')]);
title(['L= ',num2str(L), ', \mu= ',num2str(mu), ...
    ', SNR= ',num2str(SNR), ', M= ', num2str(M)]);

figure(3); plot_filter(h, h_hat(:,:,3)) % filter coeff.