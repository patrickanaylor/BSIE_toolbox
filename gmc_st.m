function [ClustMtx,ClustMem,ClustNum] = gmc_st(zr,tol)

% This function attempts to cluster elements in the "zr" matrix given the
% constriant that
% a) each member within a cluster must come from an element in each column
%    of "zr" (i.e., number of members in a cluster = size(zr,2)),
% b) that all the Euclidean distances between the members inside a cluster
%    must be less than or equal to "tol".
%
% Note that any member in a cluster can belong to other clusters
% (one-to-many relationship).
%
%	Input Parameters [size]:
%       zr       : zeros of the channel impulse responses [L-1 x M]
%       tol      : tolerence
%
%	Output Parameters:
%       ClustMtx : a matrix containing number of clusters-by-cluster
%                  members
%
%	References:
%       [1] A. W. H. Khong, X. Lin, and P. A. Naylor, 
%           "Algorithms for identifying clusters of near-common 
%           zeros in multichannel blind system identification and
%           equalization," in Proc. IEEE Intl. Conf. Acoust.,
%           Speech, Signal Processing (ICASSP), Las Vegas, USA,
%           Apr. 2008, pp. 389-392.
%
% Authors: A. Kong
%
% History: 2007-07-05 Initial version
%
% Copyright (C) Imperial College London 2009-2010

[len, totCh] = size(zr);

[pwIdxCm, pwChCm]   = compdist(zr,tol,'n');
resultTst           = chkpwChCm(pwChCm,totCh);

if totCh==2;
    [pwIdxCm, pwChCm]   = compdist(zr,tol,'n');
    ClustMtx            = pwIdxCm;
else
    if isempty(pwIdxCm)
        disp('No clusters found: terminated by zero pairwise channels.');
        ClustMtx = [];
    elseif (resultTst==1)
        [ClstrComb,MemLen]      = GrpClstrComb(pwIdxCm, pwChCm);
        minMemLen               = min(MemLen);
        minMemLenIdx            = find(MemLen==minMemLen);
        minGpIdx                = minMemLenIdx(1);
        sumchCombMtx            = ChkChComb(ClstrComb,totCh);
        
        if sumchCombMtx==totCh*(totCh-1)/2;
            ClustMtx = getClusters(ClstrComb,minGpIdx,minMemLen,totCh,tol,zr);
        else
            disp('No clusters found: terminated by missing pairwise channels.');
            ClustMtx = [];
        end
    elseif (resultTst==0);
        disp('No clusters found: terminated by missing pairwise channels.');
        ClustMtx = [];
    end
end

ClustIdx = ClustMtx;
ClustMem = getClustMem(zr,ClustIdx,totCh);
ClustNum = size(ClustIdx,1);

%--------------------------------------------------------------------------
function [ClstrComb,MemLen] = GrpClstrComb(pwIdxCm, pwChCm)
% this function tabulates cluster pair matrix according to their pairwise
% channels
absdif = abs(pwChCm-circshift(pwChCm,1));
sumdif = sum(absdif');
chIdx  = find(sumdif~=0);
for i = 1:length(chIdx);
    ClstrComb.ChGp{i}  = pwChCm(chIdx(i),:);
end
for i = 1: length(chIdx)-1;
    ClstrComb.MemGp{i} = pwIdxCm(chIdx(i):chIdx(i+1)-1,:);
    MemLen(i)          = size(ClstrComb.MemGp{i},1);
end
ClstrComb.MemGp{i+1} = pwIdxCm(chIdx(end):size(pwIdxCm,1),:);
MemLen(i+1)          = size(ClstrComb.MemGp{i+1},1);

%--------------------------------------------------------------------------
function out = ChkChComb(ClstrComb,totCh)
% this function checks that all channels are paired. If "out" is not equal
% to 0.5M(M-1), then no common zeros;
chCombMtx    = zeros(totCh,totCh);
for i = 1:length(ClstrComb.ChGp);
    chCombMtx(ClstrComb.ChGp{i}(1),ClstrComb.ChGp{i}(2)) = 1;
end
out = sum(sum(chCombMtx));

%--------------------------------------------------------------------------
function ClustMtx = getClusters(ClstrComb,minGpIdx,minMemLen,totCh,tol,zr);
% this function extracts out the clusters;
ClustMtx = [];
for i = 1:minMemLen
    seedMem             = ClstrComb.MemGp{minGpIdx}(i,:);
    seedCh              = ClstrComb.ChGp{minGpIdx};
    ClustVec            = updateClust(seedMem,seedCh,ClstrComb,totCh,tol,zr); % can be a vector or matrix (chared cluster)
    ClustMtx            = [ClustMtx; ClustVec];
end

%--------------------------------------------------------------------------
function ClustVecOut = updateClust(seedMem,seedCh,ClstrComb,totCh,tol,zr);
% this is the main fn where clusters get expanded(shared), or gets
% trimmed.
ClustVec          = zeros(1,totCh);
ClustVec(seedCh)  = seedMem;
eptyCh            = find(ClustVec==0);
noneptyCh         = find(ClustVec~=0);
ClustVecOut       = [];

while ~isempty(eptyCh);
    eptyCh  = find(ClustVec(1,:)==0);
    ChB     = eptyCh(1);
    ChA     = noneptyCh(1);
    
    memA            = ClustVec(1,ChA);
    [actMem,chIx]   = getActMem(ClstrComb,ChA,ChB);
    rww             = find(actMem(:,chIx)==memA);
    chIy            = swapIx(chIx);
    memB            = actMem(rww,chIy);
    
    if isempty(memB);
        ClustVec(1,:)=[];
    else
        ClustVecConc                 = repmat(ClustVec(1,:),size(memB,1)-1,1);
        ClustVec                     = [ClustVecConc;ClustVec];
        ClustVec(1:size(memB,1),ChB) = memB;
    end
    
    if size(ClustVec,1)>=1;
        ClustVec        = trimClustVec(ClustVec,tol,zr);
    end
    
    if size(ClustVec,1)>=1;
        numNZ            = all(ClustVec');
        NZid             = find(numNZ==1);
        ClustVecOut      = [ClustVecOut;ClustVec(NZid,:)];
        ClustVec(NZid,:) = [];
    end
    
    if size(ClustVec,1)>=1;
        eptyCh = find(ClustVec(1,:)==0);
    else
        eptyCh = [];
    end
    
end

%--------------------------------------------------------------------------
function [actMem,chIx] = getActMem(ClstrComb,ChA,ChB)
% this extracts out members of ChA and ChB
for j=1:length(ClstrComb.ChGp);
    gt1     = (ClstrComb.ChGp{j}==[ChA ChB]);
    gt2     = (ClstrComb.ChGp{j}==[ChB ChA]);
    smtst1   = sum(gt1);
    smtst2   = sum(gt2);
    if (smtst1==2)||(smtst2==2);
        ot  = j;
        chIx  = find(ClstrComb.ChGp{j}==ChA);
        break;
    end
end
actMem = ClstrComb.MemGp{ot};

%--------------------------------------------------------------------------
function out = trimClustVec(ClustVec,tol,zr)
% this function trims off ClustVec if any of the pairs are not within tol
rwNum= size(ClustVec,1);
out  = [];
for rr=1:rwNum;
    Clust       = ClustVec(rr,:);
    nzrId       = find(Clust~=0);
    nzrNum      = length(nzrId);
    pp          = NaN*ones(nzrNum,nzrNum);
    for ii=1:nzrNum;
        pp(ii,ii) = zr(Clust(nzrId(ii)),nzrId(ii));
    end
    [A,C] = compdist(pp,tol,'n');
    if size(A,1)==0.5*nzrNum*(nzrNum-1);
        out = [out;Clust];
    end
end

%--------------------------------------------------------------------------
function chIy = swapIx(chIx)
if chIx==2;chIy=1;
else chIy=2;
end

%--------------------------------------------------------------------------
function ClustMem = getClustMem(zr,ClustIdx,totCh)
ClustMem = ClustIdx.*0;

for rw= 1: size(ClustIdx,1);
    ClustMem(rw,:) = diag(zr(ClustIdx(rw,:),1:totCh))';
end

%--------------------------------------------------------------------------
function out=chkpwChCm(pwChCm,totCh)
absdif = abs(pwChCm-circshift(pwChCm,1));
sumdif = sum(absdif');
chIdx  = find(sumdif~=0);
out    = 0;
if length(chIdx)==totCh*(totCh-1)/2;
    out=1;
end