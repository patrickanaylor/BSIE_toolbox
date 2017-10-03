function phandle = plot_npm(npm_dB, fs, ns, varargin)

% Copyright (C) Imperial College London 2009-2010

if nargin < 3
    ns = 1;
end  

ax  = 1/fs : ns/fs : length(npm_dB)*ns/fs;
phandle = plot(ax, npm_dB, varargin{:}); 
grid on; 
xlabel('Time (s)'); 
ylabel('NPM (dB)'); 