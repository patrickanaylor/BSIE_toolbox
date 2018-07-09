function [emag emag_freq] = magnitude_deviation(heq, NFFT)

% Magnitude distortion measure
%
%   [emag emag_freq] = magnitude_deviation(heq, NFFT)
%
%   Input Parameters:
%       heq       :  equalized impulse response
%       NFFT      :  length discrete Fourier transform
%
%   Output Parameters:
%       emag      :  magnitude deviation
%       emag_freq :  magnitude deviation per frequency
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
        warning('BSIE:magnitude_distortion','FFT length is too small!');
    end
end

% Magnitude distortion
Heq = fft(heq,NFFT);    
Hdb = mean(10*log10(abs(Heq(1:NFFT/2+1))));
emag = sqrt(mean((10*log10(abs(Heq(1:NFFT/2+1))) - Hdb).^2));  % RMS
if nargout > 1
    emag_freq = 10*log10(abs(Heq(1:NFFT/2+1))) - Hdb;
end