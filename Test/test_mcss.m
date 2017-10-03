% Run subspace BSI algorithm: MCSS
% This script runs the MCSS algorithm for blind 
% system identification
%
% References:
%   [1] J. Allen and D. Berkley, "Image method for 
%       efficiently simulating small room acoustics",
%       JASA, 1979.
%
% Authors: E.A.P. Habets
%
% History: 2010-04-15 Initial version by EH
%
% Copyright (C) Imperial College London 2010

clear; 
close all;

%% Initialization
M = 4;      % number of channels
fs = 8e3;   % sampling frequency
L = 64;     % channel length
SNR = 50;   % signal to noise ratio (dB)
N = 2*fs;   % data length (samples)

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.2;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
h = generate_data(M, L, fs, air, N, SNR);

% theta = pi/10; % IC
% phi = pi/10; % IC
% theta = pi/10; % WC
% phi = 2*pi/3; % WC
% h(:,1) = [1 -2*cos(theta) 1]';
% h(:,2) = [1 -2*cos(theta+phi) 1]';
% h(:,3) = [1 -2*cos(theta+2*phi) 1]';

% RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 3));
% h = randn(L,M);
% h = h ./ norm(h,2);

%% Generate signals
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));
s = randi(2,1,N);   % Generate N length of either 0 or 1.
s = s - mean(s);
for i = 1:M
    x(:,i) = filter(h(:,i),1,s);
end

%% Add noise
for i = 1:M
    x(:,i) = awgn(x(:,i),SNR);
end

%% Processing Loop: run MCLS
% methods = {'Moulines'};
methods = {'Moulines','Moulines_constraint'};
npm_dB = zeros(length(methods),2);
% h_hat = zeros(L,M,length(methods));
D = 0; % Over-estimate with D coefficients
for mm = 1:length(methods)
    % Overestimate but given the correct channel order
    % h_hat(:,:,mm) = mcss(x,L+2,L,methods{mm});
    
    % Overestimate and given the incorrect channel order
    h_hat(:,:,mm) = mcss(x,L+D,L+D,methods{mm});
    
    npm_dB(mm,1) = 20*log10(npm(h, h_hat(:,:,mm)));
    npm_dB(mm,2) = 20*log10(npm_ac(h, h_hat(:,:,mm)));
end

%% Plot results
fprintf('NPM\n');
for mm = 1:length(methods)
    fprintf([methods{mm} ' : %.2f / %.2f dB\n'], npm_dB(mm,:));
end

NFFT=2^nextpow2(L+1);
figure(2);
Ptmp=20*log10(abs(fft(h(:,1),NFFT))); 
plot(Ptmp(1:NFFT/2+1)-mean(Ptmp),'g');
hold on; 
Ptmp = 20*log10(abs(fft(h_hat(:,1,1),NFFT)));
plot(Ptmp(1:NFFT/2+1)-mean(Ptmp),':r');
Ptmp = 20*log10(abs(fft(h_hat(:,1,2),NFFT)));
plot(Ptmp(1:NFFT/2+1)-mean(Ptmp),'b');
hold off; grid on;
legend('True','Moulines','Habets')
xlabel('Frequency');
ylabel('Energy (dB)');

figure;plot([[zeros(D,1) ; h(:,1)] h_hat(:,1,1) -h_hat(:,1,2)]);