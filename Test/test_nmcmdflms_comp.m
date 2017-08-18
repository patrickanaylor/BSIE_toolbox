% Run adaptive BSI algorithm: NMCMDFLMS
% This script runs the NMCMDFLMS algorithm [2] for
% blind system identification
%
% References:
%   [1] J. Allen and D. Berkley, "Image method for
%       efficiently simulating small room acoustics",
%       JASA, 1979.
%   [2] R. Ahmad, A. W. H. Khong, and P. Naylor, "Proportionate frequency
%       domain adaptive algorithms for blind channel identification," in
%       Proc. ICASSP, Toulouse, France, May 2006, vol. 5, pp. 29-32.
%
% Authors: E.A.P. Habets
%
% History: 2011-03-01 Initial version
%
% Copyright (C) Imperial College London 2011

clc;
clear;
close all;

%% Initialization
Mch = 5;            % number of channels
fs = 8e3;           % sampling frequency
L = 128;            % channel length
F = 2*L;            % frame lenght
SNR = 40;           % signal to noise ratio (in dB)
TimeSec = 10;       % Duration (seconds)
Ldata = TimeSec*fs; % data length (samples)
rho = [0.6 ; 0.6];  % step-size
lambda = [];
delta = [];

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(Mch, L, fs, air, Ldata, SNR);

%% Process
for Mk = 2
    Ms = 2^Mk;
    Ns = 2*L/Ms;
    npm_dB_NMCMDFLMS_all = [];
    npm_dB_NMCFLMS_all = [];
    
    for rho_idx = 1:size(rho,2)
        % Initiailize NMCMDFLMS
        [shMDF h_hat] = init_nmcmdflms(L, Ms, Mch, x);
        
        ns = Ns/2; % number of new samples per iteration
        B = fix(Ldata/ns); % number of input blocks of length L
        npm_dB = zeros(B,1);
        
        % Processing Loop: run NMCMDFLMS
        wbar = waitbar(0,'NMCMDFLMS');
        
        for bb = 1 : B
            waitbar(bb/B);
            
            if bb == 1
                xm = [zeros(ns,Mch); x(1:ns,:)];
            else
                xm = [xm(ns+1:end,:); x((bb-1)*ns+1:bb*ns,:)];
            end
            
            shMDF.bb = bb;
            shMDF.normVal = 1;
            [shMDF, h_hat, J_tmp] = nmcmdflms(shMDF, xm, rho(1,rho_idx), lambda, delta);
            npm_dB(bb) = 20*log10(npm(h, h_hat));
        end
        npm_dB_NMCMDFLMS_all = [npm_dB_NMCMDFLMS_all; npm_dB];
        close(wbar)
    end
end

for rho_idx = 1:size(rho,2)
    % Initiailize NMCFLMS
    [h_hat, P_k_avg] = init_nmcflms(L, F, Mch, x(1:F,:));
    
    ns = F-L;
    B = fix(Ldata/ns);  % number of input blocks of length L
    J = zeros(Ldata,1);
    npm_dB = zeros(B,1);
        
    % Processing Loop: run NMCFLMS
    wbar = waitbar(0,'NMCFLMS');
    for bb = 1 : B
        waitbar(bb/B);
        if bb == 1
            xm = [zeros(F-ns,Mch); x(1:ns,:)];
        else
            xm = [xm(ns+1:end,:); x((bb-1)*ns+1:bb*ns,:)];
        end
        [h_hat, P_k_avg, J_tmp] = nmcflms(xm, h_hat,...
            P_k_avg, rho(2,rho_idx), lambda, delta);
        J((bb-1)*ns+1:(bb)*ns) = J_tmp(F-L-ns+2:end);
        npm_dB(bb) = 20*log10(npm(h, h_hat));
    end
    npm_dB_NMCFLMS_all = [npm_dB_NMCFLMS_all; npm_dB];
    close(wbar);
end

%% Plot results
tsec = Ldata/fs;

figure(1);
hold on
xscale = (1:length(npm_dB_NMCFLMS_all))/length(npm_dB_NMCFLMS_all)*tsec;
plot(xscale,npm_dB_NMCFLMS_all,'b-');
xscale = (1:length(npm_dB_NMCMDFLMS_all))/length(npm_dB_NMCMDFLMS_all)*tsec;
plot(xscale,npm_dB_NMCMDFLMS_all,'r--');
hold off;
rho_str = num2str(rho.');
legend('NMCFLMS','NMCMDFLMS Ms=4');
grid on;
xlabel('Time (s)'); 
ylabel('NPM (dB)'); 