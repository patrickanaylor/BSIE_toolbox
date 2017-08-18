function plot_J(J, fs, varargin)

% Copyright (C) Imperial College London 2009-2010

ax  = 1/fs : 1/fs : length(J)/fs;
plot(ax, 10*log10(J),varargin{:}); grid on;
xlabel('Time (s)'); ylabel('J (dB)'); 
