function [h_hat, R_hat] = mcls(x,L,method,h_hat)

% Multichannel LS cross-relation algorithms for blind SIMO system identification
%
%	[h_hat, R_hat, J] = mcls(x,L,method,h_hat)
%
%   Input parameters:
%   	xin    : input matrix [N x M]
%       L      : channel length
%       method : {'recursive','xu','minimum','minimum_unbiased','spcc'}
%       h_hat  : updated filter coef. matrix [L x M] (only for SPCC)
%
%   Outputs parameters:
%       h_hat  : updated filter coef. matrix [L x M]
%       R_hat  : covariance matrix [M L x M L]
%
% Authors: E.A.P. Habets
%
% History: 2010-03-05 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2010

error(nargchk(2,4,nargin));

if ~exist('method','var');
    method = 'xu';
end

if strcmp(method,'recursive')
    lambda = 1;
end

M = size(x,2);
N = size(x,1);

switch lower(method);
    case 'recursive'
        R_hat = zeros(M*L,M*L);
        for nn = L : N
            xin = x(nn:-1:nn-L+1,:);
            R_tilde = -xin(:)*xin(:)';
            Rr = xin*xin';
            for jj = 1:M
                row = (jj-1)*L+1:jj*L;
                for ii = 1:M
                    col = (ii-1)*L+1:ii*L;
                    R_tilde(row,col) = R_tilde(row,col)';
                end
                R_tilde(row,row) = R_tilde(row,row) + Rr;
            end
            R_hat = lambda * R_hat + R_tilde;
        end
        
        R_hat = R_hat./((N-L)*(M-1));
        
    case 'recursive_alt' % Test by EH
        R = zeros(M*L,M*L);
        for nn = L : N
            xin = x(nn:-1:nn-L+1,:);
            R_tilde = xin(:)*xin(:)';
            for jj = 1:M
                row = (jj-1)*L+1:jj*L;
                for ii = 1:M
                    col = (ii-1)*L+1:ii*L;
                    R_tilde(row,col) = R_tilde(row,col)';
                end
            end
            R = R + R_tilde;
        end
        
        D = zeros(L*M,L*M);
        for ii = 1:M
            row = (ii-1)*L+1:ii*L;
            for nn = [1:ii-1 ii+1:M]
                idx = (nn-1)*L+1:nn*L;
                D(row,row) = D(row,row) + R(idx,idx);
            end
        end
        
        R_hat = zeros(L*M,L*M);
        for ii = 1:M
            row = (ii-1)*L+1:ii*L;
            for jj = 1:M
                col = (jj-1)*L+1:jj*L;
                if ii==jj
                    R_hat(row,col) = D(row,col);
                else
                    R_hat(row,col) = -R(row,col);
                end
            end
        end
        
    case 'xu'
        X = [];
        for ii = 1:M-1
            X_left = [];
            for ll = ii+1:M
                X_tmp = convmtx(x(:,ll),L);
                X_left = [X_left ; X_tmp(L:N-L+1,:)];
            end
            X_tmp = convmtx(x(:,ii),L);
            X_right = repblkdiag(-X_tmp(L:N-L+1,:),M-ii);
            X = [ X ; zeros((N-2*L+2)*(M-ii),L*(ii-1)) X_left X_right];
        end
        
        R_hat = X.'*X;
        R_hat = R_hat./((N-L)*(M-1));
        
    case 'convmtx' % Not OK for small N!
        X = [];
        for ii = 1:M-1         
            X_left = [];
            for ll = ii+1:M
                X_left = [X_left ; convmtx(x(:,ll),L)];
            end 
            X_right = repblkdiag(-convmtx(x(:,ii),L),M-ii);
            X = [ X ; zeros((N+L-1)*(M-ii),L*(ii-1)) X_left X_right];
        end

        R_hat = X.'*X;
        R_hat = R_hat./((N-L)*(M-1));
        
   case 'minimum'    
        all_perm = nchoosek(1:M,2);
        [~,pos] = unique(all_perm(:,2));
        
        % For nn=1
        X_1 = convmtx(x(:,1),L);
        X_2 = convmtx(x(:,2),L);
        X = [X_2(L:N-L+1,:) -X_1(L:N-L+1,:)];
        
        % For 2<=nn<=length(pos)
        for nn = 2:length(pos)
            X_tmp = convmtx(x(:,all_perm(pos(nn),2)),L);
            X_left = [zeros(N-2*L+2,L*(all_perm(pos(nn),1)-1)) X_tmp(L:N-L+1,:)];
            X_tmp = convmtx(x(:,all_perm(pos(nn),1)),L);
            X_right = -X_tmp(L:N-L+1,:);
            X = [X zeros(size(X,1),L); X_left X_right];
        end
        
        R_hat = X.'*X;
        R_hat = R_hat./((N-L)*(M-1));
        
    case 'minimum_unbiased'
        all_perm = nchoosek(1:M,2);
        all_perm = [all_perm ; M 1];
        [~,pos] = unique(all_perm(:,2));
        pos = sort(pos);
        
        % For nn=1
        X_1 = convmtx(x(:,1),L);
        X_2 = convmtx(x(:,2),L);
        X = [X_2(L:N-L+1,:) -X_1(L:N-L+1,:)];
        
        % For 2<=nn<=lenght(pos)-1
        for nn = 2:length(pos)-1
            X_tmp = convmtx(x(:,all_perm(pos(nn),2)),L);
            X_left = [zeros(N-2*L+2,L*(all_perm(pos(nn),1)-1)) X_tmp(L:N-L+1,:)];
            X_tmp = convmtx(x(:,all_perm(pos(nn),1)),L);
            X_right = -X_tmp(L:N-L+1,:);
            X = [X zeros(size(X,1),L); X_left X_right];
        end
        
        % For nn==length(pos)
        nn = length(pos);
        X_tmp = convmtx(x(:,all_perm(pos(nn),1)),L);
        X_left = [X_tmp(L:N-L+1,:) zeros(N-2*L+2,L*(all_perm(pos(nn),1)-2))];
        X_tmp = convmtx(x(:,all_perm(pos(nn),2)),L);
        X_right = -X_tmp(L:N-L+1,:);
        X = [X; X_left X_right];
                
        R_hat = X.'*X;
        R_hat = R_hat./((N-L)*(M-1));
        
    case 'spcc'        
        R = zeros(L*M,L*M);
        for nn = L : N
            xin = x(nn:-1:nn-L+1,:);
            R_tilde = xin(:)*xin(:)';
            for jj = 1:M
                row = (jj-1)*L+1:jj*L;
                for ii = 1:M
                    col = (ii-1)*L+1:ii*L;
                    R_tilde(row,col) = R_tilde(row,col)';
                end
            end
            R = R + R_tilde;
        end
        
        alpha = zeros(M,M);
        beta = zeros(M,M);
        for ii=1:M
            row = (ii-1)*L+1:ii*L;
            for jj=1:M
                col = (jj-1)*L+1:jj*L;
                alpha(ii,jj) = (h_hat(:,ii).'*R(col,row)*h_hat(:,jj))/(h_hat(:,ii).'*R(col,col)*h_hat(:,ii));
                beta(ii,jj) = alpha(ii,jj)/(h_hat(:,jj).'*R(row,row)*h_hat(:,jj));
            end
        end
        
        D = zeros(L*M,L*M);
        for ii = 1:M
            row = (ii-1)*L+1:ii*L;
            for nn = [1:ii-1 ii+1:M]
                idx = (nn-1)*L+1:nn*L;
                D(row,row) = D(row,row) + alpha(ii,nn) * beta(ii,nn) * R(idx,idx);
            end
        end
        
        R_hat = zeros(L*M,L*M);
        for ii = 1:M
            row = (ii-1)*L+1:ii*L;
            for jj = 1:M
                col = (jj-1)*L+1:jj*L;
                if ii==jj
                    R_hat(row,col) = D(row,col);
                else
                    R_hat(row,col) = -beta(ii,jj) * R(row,col);
                end
            end
        end

    otherwise
        error('Unknown method.');         
end

[~,~,V] = svd(R_hat);
h_hat = V(:,end) / norm(V(:,end));
h_hat = reshape(h_hat,L,M);

% [V,DD] = eig(R_hat);
% [~,pos] = min(diag(DD));
% h_hat = reshape(V(:,pos),L,M);
% h_hat = h_hat./norm(h_hat);

% figure(1);imagesc(20*log10(abs(X)))

function B = repblkdiag(A,num)
B = [];
for l = 1:num
    B = blkdiag(B,A);
end