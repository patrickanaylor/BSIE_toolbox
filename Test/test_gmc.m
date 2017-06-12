% Runs the GMC_ST algorithms to find
% the number of near-common zeros
%
% References:
%   [1] J. Allen and D. Berkley, "Image method for
%       efficiently simulating small room acoustics",
%       JASA, 1979.
%
% Authors: E.A.P. Habets
%
% History: 2009-07-14 Initial version
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: test_gmc.m 425 2011-08-12 09:15:01Z mrt102 $

clc;
clear;
close all;

%% Initialization
M = 2;       % number of channels
fs = 8e3;    % sampling frequency
L = 512;     % channel length
tol = 0.005; % tolerance

% Generate AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
h = generate_data(M, L, fs, air);

% Find zeros of the acoustic impulse responses
for ii=1:M;
    zr(:, ii) = roots(h(:, ii));
end

% Find number of near-common zeros using GMC-DC
% display(sprintf('\nRunning GMC-DC...'))
% ZeroMtx = gmc_dc(zr, tol);
% display(sprintf('GMC-DC: %d\n', size(ZeroMtx, 1)));
% clear ZeroMtx;

% Find clusters of near-common zeros using GMC-ST
display(sprintf('\nRunning GMC-ST...'))
[ClustMtx, ClustMem, ClustNum]  = gmc_st(zr, tol);
display(sprintf('GMC-ST: %d',ClustNum));