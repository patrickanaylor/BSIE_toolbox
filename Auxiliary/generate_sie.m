function ie = generate_sie(h, npm, mode, T60, fs)

% This function generates system identification errors
%
%   [ie]  = generate_sie(h, npm, mode, T60, fs)
%
%	Input Parameters:
%       h      : system impulse responses
%       npm    : desired NPM level in dB
%       mode   : the shape of errors. Three options are 
%                'wgn'  white Gaussian noise like errors;
%                'prop' WGN but proportional to h;
%                'damp' damping errors, damping factor same as h
%       T60    : the reverberation time of h, only needed for
%                generating damping errors
%       fs     : sampling frequency, only needed for generating damping
%                errors (optional: default 8000)
%
%	Output Parameters:
%       ie     : the representation of identification errors
%
% Author : W. Zhang
%
% History: 2009-07-13 Initial version
%
% Copyright (C) Imperial College London 2009-2010

[L M] = size(h);
re = zeros(L, M);
ie = zeros(L, M);

if nargin == 3 && strcmp(mode, 'damp')
    error('T60 is required to generate damping errors');
else
    for m = 1:M
        switch mode
            case 'wgn'
                re(:,m) = randn(L,1);
            case 'prop'
                re(:,m) = randn(L,1) .* h(:,m);
            case 'damp'
                if nargin == 4
                    fs = 8000;
                end
                alpha = 3*log(10)/T60/fs; 
                re(:,m) = randn(L,1) .* exp(-alpha*(1:L))';
        end
        ie(:,m) = re(:,m)-(re(:,m)'*h(:,m))/(h(:,m)'*h(:,m))*h(:,m);
        ie(:,m) = ie(:,m)/norm(ie(:,m))*norm(h(:,m))*10^(npm/20)/sqrt(1-10^(2*npm/20));
    end
end