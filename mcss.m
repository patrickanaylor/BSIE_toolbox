function [h_hat, R_hat] = mcss(x,L,Lc,method)

% Multichannel subspace algorithms for blind SIMO system identification
%
%	[h_hat, R_hat, J] = mcss(x,L,Lc,method)
%
%   Input parameters:
%   	xin    : input matrix [N x M]
%       L      : overestimated channel length (Note: Moulines uses N)
%       Lc     : channel length (Note: Moulines uses M=Lc-1)
%       method : {'Moulines'}
%
%   Outputs parameters:
%       h_hat  : updated filter coef. matrix [L x M]
%       R_hat  : covariance matrix [M L x M L]
%
% Authors: E.A.P. Habets
%
% History: 2010-04-15 Initial version by E.A.P. Habets
%
% Copyright (C) Imperial College London 2010

narginchk(2,4);

if ~exist('method','var');
    method = 'moulines';
end

M = size(x,2);
N = size(x,1);

switch lower(method);
    case 'moulines'
        R_hat = zeros(M*L,M*L);
        for nn = L : N
            xin = x(nn:-1:nn-L+1,:);
            Rr = xin(:)*xin(:)';
            R_hat = R_hat + Rr;
        end
        R_hat = R_hat./((N-L)*M);
        
        [~,~,V] = svd(R_hat);

        V_null = V(:,L+Lc:end);

        for k = 1:size(V_null,2) % M*L-(Lc-1)-L
            Gi = MC_FilteringMatrix(V_null(:,k),L,Lc-1);    
            if k ==1
                Q = Gi*Gi';
            else
                Q = Q + Gi*Gi';
            end
        end
        [~,~,V] = svd(Q);

        h_hat = reshape(V(:,end),Lc,M);
        h_hat = h_hat./norm(h_hat);     
        
    case 'moulines_constraint' % Under construction... (EH)      
        R_hat = zeros(M*L,M*L);
        for nn = L : N
            xin = x(nn:-1:nn-L+1,:);
            Rr = xin(:)*xin(:)';
            R_hat = R_hat + Rr;
        end
        R_hat = R_hat./((N-L)*M);
        
        [~,~,V] = svd(R_hat);

        V_null = V(:,L+Lc:end);

        for k = 1:size(V_null,2) % M*L-(Lc-1)-L
            Gi = MC_FilteringMatrix(V_null(:,k),L,Lc-1);    
            if k ==1
                Q = Gi*Gi';
            else
                Q = Q + Gi*Gi';
            end
        end
        [~,S,V] = svd(Q);
        
        K = diff(log(diag(S)));
        [~,pos_min] = min(K);
        % display([num2str(pos_min) '/' num2str(length(diag(S)))]);
        
        % Maximize passband engery!
        V_prime = V(:,pos_min+1:end);       
        F = zeros(Lc,Lc);
        BP = 4;
        for k = BP:ceil(Lc/2+1)-BP
            f = exp(-1i*2*pi*(0:Lc-1)*k/Lc).';
            F = F + f * f';
        end
        B = repblkdiag(F,M);       
        [~,~,V] = svd(V_prime' * B * V_prime); 
        beta = V(1,:).';
        
        if isreal(x)
            h_hat = real(reshape(V_prime*beta,Lc,M));
        else
            h_hat = reshape(V_prime*beta,Lc,M);
        end
        h_hat = h_hat./norm(h_hat);     
        
    otherwise
        error('Unknown method.');        
end

function B = repblkdiag(A,num)
B = [];
for l = 1:num
    B = blkdiag(B,A);
end

function [V_tot] = MC_FilteringMatrix(V, N, M)
[LN,~] = size(V);
L = LN/N; % Number of channels
V_tot = [];
for l = 0:L-1
    Vl = V(l*N+1:(l+1)*N).';
    Vl_extent = [Vl , zeros(1,M)];
    Vl_mat = toeplitz([Vl_extent(1),zeros(1,M)],Vl_extent);
    V_tot = [V_tot ; Vl_mat]; %#ok<AGROW>
end