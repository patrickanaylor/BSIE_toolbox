function phandle = plot_npm(npm_dB, fs, ns, varargin)

% Copyright (C) Imperial College London 2009-2010
% Version: $Id: plot_npm.m 425 2011-08-12 09:15:01Z mrt102 $

if nargin < 3
    ns = 1;
end  

ax  = 1/fs : ns/fs : length(npm_dB)*ns/fs;
phandle = plot(ax, npm_dB,varargin{:}); 
grid on; 
xlabel('Time (s)'); 
ylabel('NPM (dB)'); 