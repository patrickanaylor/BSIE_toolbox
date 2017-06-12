function [npm_val] = npm_ac(h, hhat)

% Normalized Projection Misalignment with Alignment Correction
%
%   [npm_val] = npm_ac(h, hhat)
%
%	Input Parameters [size]:
%       h       : true impulse responses [L x M]
%       hhat    : estimated impulse responses [L x M]
%
%	Output Parameter:
%       npm_val : Normalize Projection Misalignment
%
%	References:
%       [1] D. R. Morgan, J. Benesty and M. M. Sondhi, "On the evaluation of
%           estimated impulse responses," IEEE Signal Processing Lett., Vol. 5, No.
%           7, pp. 174-176 Jul 1998.
%
% Authors: E.A.P. Habets 
%
% History: 2010-04-16 - added alignment correction to npm.m (EH)
%
% Copyright (C) Imperial College London 2010
% Version: $Id: npm_ac.m 425 2011-08-12 09:15:01Z mrt102 $

if size(hhat,1) <= size(h,1)
    hv = reshape(h(1:size(hhat,1),:),[],1);
    hhatv = hhat(:);
else
    % hv = h(:);
    % hhatv = reshape(hhat(1:size(h,1),:),[],1);
    h = [h ; zeros(size(hhat,1)-size(h,1),size(h,2))];
    hv = h(:);
    hhatv = hhat(:);
end

hhatv_sc = sign((hv.'*hhatv)/(hhatv.'*hhatv))*hhat(:);
hhat_sc = reshape(hhatv_sc,size(hhat,1),size(hhat,2));

% Correct possible phase shift!
L = max(size(hhat,1),size(h,1));
NFFT=2^(nextpow2(L)+1);
H = fft(h,NFFT);
Hhat = fft(hhat_sc,NFFT);
binFreq = ((mod(((0:NFFT-1)+floor(NFFT/2)), NFFT)-floor(NFFT/2))/NFFT).';
u = Hhat .* conj(H);
if mod(NFFT, 2) == 0
    u(length(u)/2+1,:) = 0;
end
Xcor = abs(ifft(u));
intDelay = find(Xcor==max(max(Xcor)));
intDelay = intDelay(1)-1;
rotN = repmat(exp(1i*2*pi*intDelay .* binFreq),1,size(h,2));
uDelayPhase = repmat(-2*pi*binFreq,1,size(h,2));
u = u .* rotN;
weight = abs(u); 
uDelayPhase = uDelayPhase .* weight;
weighted_angle_u = angle(u).* weight;
% Note: fracDelay minimizes norm(uDelayPhase(:)*fracDelay - weighted_angle_u(:),2)
fracDelay = pinv(uDelayPhase(:))*weighted_angle_u(:); 
Hhat = Hhat .* repmat(exp(1i*2*pi*(intDelay+fracDelay).*binFreq),1,size(h,2));
hhat_new = real(ifft(Hhat));
hhat = hhat_new(1:size(hhat,1),:);
hhatv = hhat(:);

% display(['Delay:' num2str(intDelay+fracDelay)]);
% H_angle = unwrap(angle(H(1:NFFT/2+1,:)));
% Hhat_angle = unwrap(angle(Hhat(1:NFFT/2+1,:)));
% Hhat_angle_new = unwrap(angle(Hhat(1:NFFT/2+1,:)));
% figure(10);plot([(H_angle(1:NFFT/2+1,2)) (Hhat_angle(1:NFFT/2+1,2)) (Hhat_angle_new(1:NFFT/2+1,2))]);

epsilon = hv-((hv.'*hhatv)/(hhatv.'*hhatv))*hhatv;
npm_val = norm(epsilon)/norm(hv);

%figure(100);
%plot([hv ((hv.'*hhatv)/(hhatv.'*hhatv))*hhatv]);