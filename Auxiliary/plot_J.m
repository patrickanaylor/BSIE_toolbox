function plot_J(J, fs, varargin)

% Copyright (C) Imperial College London 2009-2010
% Version: $Id: plot_J.m 425 2011-08-12 09:15:01Z mrt102 $

ax  = 1/fs : 1/fs : length(J)/fs;
plot(ax, 10*log10(J),varargin{:}); grid on;
xlabel('Time (s)'); ylabel('J (dB)'); 
