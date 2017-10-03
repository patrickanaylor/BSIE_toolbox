function plot_time(x,fs,varargin)
t = 0:1/fs:(length(x)-1)/fs;
if nargin == 2
    plot(t,x);
else
    plot(t,x,varargin{:});
end
xlabel('Time (s)');
ylabel('Amplitude');
axis tight;