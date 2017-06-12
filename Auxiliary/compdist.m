function [IdxCm, chCm, dist]= compdist(zr,tol,opdist)

% This function computes the distances between all the zeros for any
% arbitrary number of channels
%
%   Input   zr:     a N-by-M matrix where N= zeros of a transfer fn and
%                   M= channels
%
%           tol:    specify the tolerance of the common zero.
%           opdist: output dist(y or n)
%
%   Output  dist: a 4D matrix containg the distances where
%                 indices 1 and 2 contain all possible channels
%                 indices 3 and 4 contain dist. of all possible zeros
%                 of the tr. fn.
%
%           IdxCm: a matrix containing the index of the zr mtx which are
%                  common (total rows= toal common zeros, columns=2);
%
%           chCm:  a matrix containing the channels of the zr mtx which are
%                  common (total rows= toal common zeros, columns=2);
%
%
%   E.g.   dist(1,2,3,4) will return the distance between 3rd and 4th zero
%          of channel 1 and 2 respectively
%
%   Method: This uses the relationship: (A-B)^2= A^2-2AB+B^2;
%
%   Future versions:
%       - to increase efficiency, reduce the number of loops
%         from ch^2 to ch*(ch-1)/2
%
% History: 14th Feb 07
%          16th Feb 07- Added common zeros extraction
%          23th Feb 07- adjusted tolerance sensitivity
%          21st May 07- changed order of comparison using struct prog
%          23rd May 07- removed eps tolerence @ 80
%          09th Jul 07- include dist output option
%
% Copyright (C) Imperial College London 2009-2010
% Version: $Id: compdist.m 425 2011-08-12 09:15:01Z mrt102 $

[len,TotCh]    = size(zr);

if lower(opdist)=='n';
    dist    = zeros(len,len);
    ch1Bg   = zeros(len,len);
    ch2Bg   = zeros(len,len);
    
    
    IdxCm = [];
    chCm  = [];
    
    for chA= 1:TotCh;
        for chB= 1:TotCh;
            
            ch1 = zr(:,chA);
            ch2 = zr(:,chB);
            
            ch1Bg   = ch1*ones(1,len);
            ch2Bg   = ones(len,1)*ch2.';
            
            dist =...
                sqrt(abs(ch1Bg.^2-(2*ch1*ch2.')+ch2Bg.^2));
            
            clear ch1Bg;
            clear ch2Bg;
            if chA<chB;
                [IdxCmTp,chCmTp] = ...
                    xtrctCm(chA,chB,dist,IdxCm,chCm,tol);
                IdxCm = [IdxCm; IdxCmTp];
                chCm  = [chCm; chCmTp];
            end
            clear IdxCmTp;
            clear chCmTp;
            clear dist;
        end
    end
    dist = [];
    
    
else
    dist    = zeros(TotCh.^2,TotCh.^2,len,len);
    IdxCm = [];
    chCm  = [];
    
    for chA= 1:TotCh;
        for chB= 1:TotCh;
            
            ch1 = zr(:,chA);
            ch2 = zr(:,chB);
            
            ch1Bg   = ch1*ones(1,len);
            ch2Bg   = ones(len,1)*ch2.';
            
            dist(chA,chB,:,:) =...
                sqrt(abs(ch1Bg.^2-(2*ch1*ch2.')+ch2Bg.^2));
            
            clear ch1Bg;
            clear ch2Bg;
            
            if chA<chB;
                [IdxCmTp,chCmTp] = ...
                    xtrctCm(chA,chB,squeeze(dist(chA,chB,:,:)),IdxCm,chCm,tol);
                IdxCm = [IdxCm; IdxCmTp];
                chCm  = [chCm; chCmTp];
            end
            
            clear IdxCmTp;
            clear chCmTp;
        end
    end
end

% find index of common zeros and their respective channels
function [IdxCm,chCm]= xtrctCm(chA,chB,data,IdxCm,chCm,tol)
IdxCm   = [];
chCm    = [];
indx2   = [];
IdxCmI  = [];

for i= 1:size(data,1)
    indx2  = [indx2 struct('lesTolIdx',find(data(i,:)<=tol))];
end

for i=1:length(indx2);
    indx2tp = indx2(i).lesTolIdx;
    indx1tp = repmat(i,length(indx2tp),1);
    IdxCmTp  = [IdxCmI; indx1tp indx2tp'];
    if ~isempty(IdxCmTp);
        len         = size(indx1tp,1);
        chCmTp      = [chA chB];
        chCmTp      = repmat(chCmTp,len,1);
        IdxCm       = [IdxCm; IdxCmTp];
        chCm        = [chCm; chCmTp];
    end
end