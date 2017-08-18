function plot_filter(h, h_hat)

% Copyright (C) Imperial College London 2009-2010

h_norm = h(:, 1)./max(abs(h(:, 1)));
h_hat_norm = h_hat(:, 1)./max(abs(h_hat(:, 1)));

stem(h_norm,'b'); 
hold on; 
stem(h_hat_norm, 'r'); 
grid on; 
hold off;
xlabel('Time (samples)'); ylabel('Amplitude');
legend('h','h_{hat}');
axis tight;