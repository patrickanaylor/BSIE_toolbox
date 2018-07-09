function [epha epha_freq] = phase_deviation(heq, NFFT)

% Phase distortion measure
%
%   [epha epha_freq] = phase_deviation(heq, NFFT)
%
%   Input Parameters:
%       heq       :  equalized impulse response
%       NFFT      :  length discrete Fourier transform
%
%   Output Parameters:
%       epha      :  phase deviation
%       epha_freq :  phase deviation per frequency
%
% Authors: N.D. Gaubitch, E.A.P. Habets
%
% History: 2006-09-15 Initial Version by  N.D. Gaubitch
%          2007-03-12 Minor changes by N.D. Gaubitch
%          2009-07-10 NFFT option by E.A.P. Habets
%
% Copyright (C) Imperial College London 2006-2010

% Initialization
L = length(heq);
if nargin == 1
    NFFT = 2^nextpow2(2*(L-1));
else
    if NFFT < 2*(L-1)
        warning('BSIE:phase_distortion','FFT length is too small!');
    end
end

% Phase distortion
Heq = fft(heq,NFFT);    
xdat = (0:NFFT/2).';
Hang = unwrap(angle(Heq(1:NFFT/2+1)));
gamma = polyfit(xdat,Hang,1);
Hlin = gamma(1)*xdat + gamma(2);
epha = sqrt(mean((Hang-Hlin).^2));  % RMS   

if nargout > 1
    epha_freq = Hang-Hlin;
end