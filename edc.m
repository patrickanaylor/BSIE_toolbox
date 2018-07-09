function [decay, varargout] = edc(h, fs, varargin)
%   Calculate the energy decay curve and derived measures for an 
%   impulse response according to [1]
%
%   [decay,...] = edc(h,fs,...)
%
%   Inputs:
%       h       Lx1 vector containing the room impulse response
%       fs      [Optional] Sampling frequency (Hz) - required for meas
%       meas    [Optional] Any number of strings containing
%                   'Txx': Reverberation time (ms)
%                   'Cxx': Clarity index (early-to-late ratio) (dB)
%                   'Dxx': Deutlichkeit (early-to-total sound energy) (dB)
%               where xx specifies the argument in dB for 'Txx' and in ms
%               for 'Cxx' and 'Dxx'. Meas can also contain cell arrays.
%
%   Outputs:
%       decay   Lx1 normalized energy decay curve in dB
%       meas    Measures in order of input argument
%
%   Examples:
%       edc(h);
%       decay = edc(h);
%       [decay T60 C80 D50] = edc(h,8000,'T60','C80','D50');
%       [decay manymeas T60]= edc(h,8000,{'C50, 'D50'},'T60');
%
%   Notes:
%       If no outputs are specified, EDC is plotted.
%
%   To do:
%       Txx estimation could be better calculated with a linear fit.
%
%   Authors:    M. R. P. Thomas, N. D. Gaubitch
%
%   History:    26 Feb 2004 Initial version by N. D. Gaubitch
%               21 Sep 2010 Measures added by M. R. P. Thomas
%
%   References:       
%       [1] M.R. Schroeder, "New Method of Measuring Reverberation Time," 
%       J. Acoust. Soc. Am., vol. 37 pp 409-412, 1965.
%
% Copyright (C) Imperial College London 2004-2010

% Calculate EDC
h = flipud(h(:).^2);
decaylin = flipud(cumsum(h));
decaylin=decaylin/decaylin(1);
decay = 10*log10(decaylin);

% Call meas function
for ii=1:length(varargin)
    if(iscell(varargin{ii}))
        cellsize=size(varargin{ii});
        varargout{ii}=cellfun(@procmeas,repmat({decaylin},cellsize),repmat({fs},cellsize),varargin{ii});
    else
        varargout{ii}=procmeas(decaylin,fs,varargin{ii});
    end
end

% Plot if no output specified
if(nargout==0)
    if(nargin>1)
        t=1000*(0:1/fs:(length(decay)-1)/fs);
        plot(t,decay);
        xlabel('Time (ms)');
    else
        plot(decay);
        xlabel('Time (samples)');
    end
    ylabel('EDC (dB)');
    grid on;
end

function meas = procmeas(decaylin,fs,meastype)
measarg = str2double(meastype(2:end));
switch meastype(1)
    case 'T'
        [~,meas] = min(abs(10*log10(decaylin)+measarg));
        meas=1000*meas/fs;
    case 'C'
        ne = round(fs*measarg/1000);
        meas=10*log10((decaylin(1)-decaylin(ne+2))/decaylin(ne+2));
    case 'D'
        ne = round(fs*measarg/1000);
        meas=10*log10((decaylin(1)-decaylin(ne+2))/decaylin(1));    
end
