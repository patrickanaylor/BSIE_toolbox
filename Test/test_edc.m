% Runs energy decay curve (EDC) and obtains derived measures.
%
% Authors: M. R. P. Thomas
%
% History: 2010-11-29 Initial version
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: test_edc.m 425 2011-08-12 09:15:01Z mrt102 $

clc
clear
close all

%% Initialization
M = 1;     % number of channels
fs = 8e3;  % sampling frequency
L = 2048;   % channel length

% Generate AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [10; 10; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 3];  % source location
air.cen_pos = [5; 4; 1.6];  % centre pos. of the array
h = generate_data(M, L, fs, air);

mymeas = {'T30','C50','D50'};
[decay, meas] = edc(h,fs,mymeas);

t=1000*(0:1/fs:(length(decay)-1)/fs);
plot(t,decay);
xlabel('Time (ms)');
ylabel('EDC (dB)');
grid on;

for ii=1:length(mymeas)
    fprintf('%s: %f\n%',mymeas{ii},meas(ii));
end