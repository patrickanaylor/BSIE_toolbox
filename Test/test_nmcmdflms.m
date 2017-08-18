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
TimeSec = 10;       % duration (seconds)
Ldata = TimeSec*fs; % data length (samples)
rho = 0.6;
Mks = 0:2;          % determines the number of partitions to consider
lstyle = {'-','--',':'};

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x] = generate_data(Mch, L, fs, air, Ldata, SNR);

%% Process
for Mk_idx = 1:length(Mks)
    Ms = 2^Mks(Mk_idx);
    Ns = 2*L/Ms;
    npm_dB_all = [];
    
    for rho_idx = 1:length(rho)
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

            [shMDF, h_hat] = nmcmdflms(shMDF, xm, rho(rho_idx));
            npm_dB(bb) = 20*log10(npm(h, h_hat));
        end
        npm_dB_all = [npm_dB_all, npm_dB];
        close(wbar)
    end
    
    NN = length(npm_dB);
    xscale = 1:NN;
    tsec = Ldata/fs;
    xscale = xscale/NN*tsec;
    
    figure(1);
    hold on;
    plot(xscale,npm_dB_all,lstyle{Mk_idx});
end

%% Markup plot
leg = [];
for Mk_idx = 1:length(Mks)
    rho_str = num2str(rho.');
    for kk = 1:size(rho_str,1)
        leg = [leg ; ['NMCMDFLMS \rho = ', rho_str(kk,:), ', Ms = ', num2str(2^Mks(Mk_idx))]];
    end
end
grid on;
axis([0 TimeSec -25 0])
legend(leg)
hold off;
xlabel('Time (s)'); 
ylabel('NPM (dB)'); 