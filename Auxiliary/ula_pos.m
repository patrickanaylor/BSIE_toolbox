function micPos = ula_pos(micSpc,micNum,CenPos)

% This function positions a ULA such that its centroid is positioned
% at CenPos.
%
%   micPos = ula_pos(micSpc,micHght,micNum,CenPos)
%
%   Input Parameters:
%       micSpc  : microphone spacing
%       micNum  : number of microphones
%       CenPos  : center position of the array
%
%   Output parameters:
%       micPos  : a 3-by-micNum vector containing mic positions 
%
% Authors: A. Kong
%
% History: 2007 Inital version by A. Kong
%
% Copyright (C) Imperial College London 2007-2010

linArray = 0:micSpc:micSpc*(micNum-1);  % generate linear array
arrayCen = linArray(end)/2;  % find array centroid
linArray = linArray-arrayCen;  % center the array at origin
micPos(1,:) = linArray+repmat(CenPos(1),1,micNum);
micPos(2,:) = repmat(CenPos(2),1,micNum);
micPos(3,:) = CenPos(3);
