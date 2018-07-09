function [h, x, tdoa, s, src_pos, mic_pos, z] = generate_data(M, L, fs, air, N, SNR, s, v)

% This function computes the room impulse responses
% and generates WGN data with additive noise.
%
%   [h, x, tau] = generate_data(M, L, fs, air, N, SNR)
%
%	Input Parameters:
%       M      : number of channels
%       L      : filter length
%       fs     : sample frequency
%       air    : struct with acoustic impulse response parameters
%       N      : data length
%       SNR    : average signal to noise ratio across
%                sensor signals (optional: default 30 dB)
%       s      : [Optional] source signal. If absent, s generated without 
%                setting random seed. If scalar, WGN generated with seed s, 
%                else used as the source directly. Empty matrix equivalent 
%                to absent argument.               
%       v      : [Optional] noise signal. If absent, v generated without
%                setting random seed. If scalar, WGN generated with seed v, 
%                else if NxM matrix then used as the noise directly.               
%
%	Output Parameters:
%       h      : acoustic impulse responses
%       x      : noisy and reverberant sensor signals
%       tdoa   : TDOAs in (fractional) samples measured w.r.t.
%                the microphone with the shortest source-microphone
%                distance [M x 1]
%       s      : Clean signal
%       src_pos: Position of source [3x1];
%       mic_pos: Position of receivers [3xM];
%       z:     : Clean reverberant signal
%
% Author : E.A.P. Habets
%
% History: 2009-07-11 Initial version by E.A.P. Habets
%          2010-10-16 Automatically selects the microphone that is
%                     closest to the source.
%          2010-11-29 MRPT: Noise and source seeds set by argument.
%          2011-03-01 RIRs are now computed correctly.
%
% Copyright (C) Imperial College London 2009-2011

narginchk(4,8);

if ~isfield(air,'rcv_pos')
    % Determine source and sensor positions
    mic_pos = ula_pos(air.mic_spacing, M, air.cen_pos);

    [src_pos_x, src_pos_y, src_pos_z] = sph2cart(air.src_pos(1), air.src_pos(2), air.src_pos(3));
    src_pos = [mic_pos(1,1)+src_pos_x; mic_pos(2,1)+src_pos_y; mic_pos(3,1)+src_pos_z]; % speaker position
else
    mic_pos = air.rcv_pos.';
    src_pos = air.src_pos.';
end

% Compute the time of arrivals of the direct sound.
toa = zeros(M,1);
for m = 1:M
    toa(m) = norm(mic_pos(:,m)-src_pos,2)*fs/air.c;
end

% Simulate RIR using the image method [1]
L_prime = L + ceil(min(toa));
if air.T60 == 0
    h_unc = rir_generator(air.c, fs, mic_pos.', src_pos.', air.room_dim.', 1, L_prime,'omnidirectional',0,3,[],1).';
else
    h_unc = rir_generator(air.c, fs, mic_pos.', src_pos.', air.room_dim.', air.T60, L_prime,'omnidirectional',-1,3,[],1).';
end

% Align impulse response to remove direct-path propagation
% to sensor with the shortest source-microphone distance. 
% This is for the convenience of computing the NPM value.
h = rdp(h_unc,min(toa));

if nargout > 2
	tdoa = toa - min(toa);
end

% Truncate room impulse response
h = h(1:L,:);

if nargin > 4 && nargout > 1
    if nargin < 6
        SNR = 30;
    end
    
    % Generate source signal
    if ~exist('s','var')
        s = randn(N, 1);
    else
        if(length(s)==1)
            RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', s));
            s = randn(N, 1);
        elseif(isempty(s))
            s=randn(N,1);
        else
            N=length(s);
        end
    end
    
    % Calculate noisy sensor signals with correct SNR
    z = fftfilt(h, s);
    if isinf(SNR)
        x = z;
    else
        if (~exist('v','var'))
            v = randn(N,M);
        elseif(length(v)==1)
            RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', v));
            v = randn(N, M);
        end
        v = sqrt(var(s)*norm(h(:),2).^2 / (10^(SNR/10)*M*mean(var(v)))).*v;
        x = z + v;
    end
end

% DEBUG
% figure(100);
% room_plot_all(src_pos',mic_pos',air.room_dim');