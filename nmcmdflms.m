function [shMDF, H_hat, e_iijj] = nmcmdflms(shMDF, xm, rho, lambda, delta)

% The blind normalized multichannel MDF LMS (NMCMDFLMS) [1]
%
%   [shMDF H_hat J] = mcmdflms(shMDF, xm, rho, lambda, delta)
%
%   Input Parameters [size]:
%       shMDF    : MDF structure
%       xm       : input matrix [F x M]
%       H_hat    : current filter coef. matrix [L x M]
%       rho      : step size (optional: default rho=0.2)
%       lambda   : exponential forgetting factor (0 < lambda < 1) (optional)
%       delta    : regularization (optional)
%
%   Output Parameters:
%       shMDF    : MDF structure
%       H_hat    : updated filter coef. matrix [L x M]
%       J        : cross-relation error [Ns x M x M]
%
%   shMDF structure:
%       H_hat    : MDF representation of H_hat (MC)
%       Xm       : MDF representation of x - with overlapping (MC)
%       Ms       : number of MDF segments (blocks).
%       Mch      : number of channels
%       Nm       : segment length
%
%   Remarks:
%       This algorithm is suitable for estimating IRs generated using the
%       method of images.
%
%       For simulated IRs, leading zeros should be removed before employing
%       this algorithm, i.e., identify the first direct path and remove all
%       values before it and apply same process to other channels.
%
%   References:
%   [1] R. Ahmad, A. W. H. Khong, and P. Naylor, "Proportionate frequency
%       domain adaptive algorithms for blind channel identification," in
%       Proc. ICASSP, Toulouse, France, May 2006, vol. 5, pp. 29-32.
%
% Authors: B. Castro, S. Gannot and E.A.P. Habets
%
% History: 2011-03-01 Initial version
%
% Copyright (C) Bar-Ilan University 2011

% Initialization
F = size(xm,1);
U = zeros(F,shMDF.Mch,shMDF.Ms);

H_hat_10 = shMDF.H_hat;

if nargin < 3 || isempty(rho)
    rho = 0.2;
end
if nargin < 4 || isempty(lambda)
    lambda  = (1-1/(3*shMDF.Ns*shMDF.Ms/2))^(shMDF.Ns/2);
end 

Xm = fft(xm);
shMDF.Xm(:,:,2:end) = shMDF.Xm(:,:,1:end-1);
shMDF.Xm(:,:,1) = Xm;

% Calculate cross-relation error
e_iijj = zeros(shMDF.Ns,shMDF.Mch,shMDF.Mch);

% Calculate current error
for ii = 1:shMDF.Mch  % for each channel
    for jj = 1:shMDF.Mch  % CR with other channel
        for kk = 1:shMDF.Ms
            Eij =  shMDF.Xm(:,ii,kk).*H_hat_10(:,jj,kk)- ...
                shMDF.Xm(:,jj,kk).*H_hat_10(:,ii,kk);
            e_iijj(:,ii,jj) = e_iijj(:,ii,jj) - Eij;
        end
        eij = real(ifft(e_iijj(:,ii,jj)));
        eij(1:shMDF.Ns/2,:) = 0;
        e_iijj(:,ii,jj) = fft(eij);
    end
end

% Update normalization diagonal matrix
P_x = conj(Xm(:,:,1)).*Xm(:,:,1);
P_k = sum(P_x,2)*ones(1,shMDF.Mch)-P_x;
shMDF.P_k_avg = lambda*shMDF.P_k_avg + (1-lambda)*P_k;

if nargin < 5 || isempty(delta)
   delta = (shMDF.Mch-1)*mean(mean(P_x));
end

% Update filter
for ii = 1:shMDF.Mch  % for each channel
    for kk = 1:shMDF.Ms
        grad_temp = zeros(shMDF.Ns,1);
        for jj = 1:shMDF.Mch
            ttt =  conj(shMDF.Xm(:,jj,kk)) .* (e_iijj(:,ii,jj));
            grad_temp = grad_temp + ttt;
        end
        grad_temp = (1./(shMDF.P_k_avg(:,ii) + delta)).*grad_temp;
        grad_temp = real(ifft(grad_temp));
        TU = grad_temp(1:end/2,:,:);
        U(:,ii, kk) = fft([TU; zeros(size(TU))]);
    end
end

H_hat_1x = H_hat_10 - rho * U;

shMDF.H_hat = H_hat_1x;

% Normalization in time domain - to avoid trivial solution
% Firstly, calculate time domain filter from the MDF representation
% Secondly, transform the normalized time domain filter back to the MDF

shMDF_tmp = shMDF;  % temporary structure
[shMDF_tmp y_out normVal] = MDFconvCh(shMDF_tmp);
H_hat = y_out/normVal;
shMDF.normVal = normVal;
H1 = zeros(shMDF_tmp.Ns , shMDF_tmp.Mch , shMDF_tmp.Ms);
for nn = 1:shMDF_tmp.Mch
    hh1 = reshape(H_hat(:,nn),shMDF_tmp.Ns/2,shMDF_tmp.Ms);
    H = fft(hh1,shMDF_tmp.Ns);
    H1(:,nn,:) = H;
end
shMDF.H_hat = H1;

function [sMDF y_out normVal] = MDFconvCh(sMDF)
Ns = sMDF.Ns;
Ms = sMDF.Ms;
y_out = [];
x1 = [zeros(sMDF.Ns/2,1);1;zeros(Ns*Ms/2,1)];
sMDF.Xm = 0*sMDF.Xm;
B = floor(length(x1)/(Ns/2));
Ls = Ns/2;
for nn = 0:B-2
    x_s = x1(nn*Ls +1:nn*Ls + Ns);
    x_s = repmat(x_s,1,sMDF.Mch);
    
    [sMDF y_tmp] = MDFconv(sMDF,x_s);
    y_out = [y_out; y_tmp];
end
normVal = norm(y_out(:));

function [sMDF y] = MDFconv(sMDF,xnew)
sMDF.Xm(:,:,2:end) = sMDF.Xm(:,:,1:end-1);
sMDF.Xm(:,:,1) = fft(xnew);
Y = zeros(sMDF.Ns,sMDF.Mch);
for kk = 1:sMDF.Ms
    Y = Y + sMDF.Xm(:,:,kk).*sMDF.H_hat(:,:,kk);
end
yy = real(ifft(Y));
y = yy(end/2+1:end,:);