function maskPos = makeGrid(paramStruct)

% function maskPos = makeGrid(paramStruct)
%
% This function generates  a NX2 matrix of x,y positions for a grid in
% arena spatial coordinates (bottom-left pixel is <1,1>). input checking is done
% already in the higher function (generateMaskPositions). The function does
% not transforms spatial to controller coordinates (done in
% getArenaMaskTransform)
%
% INPUT
% paramStruct - structure. should contain the fields:
%   .gridSize       (1X2 vector) 
%   .startPos       (bottom-left corner of the grid - 1X2 vector) 
%   .spacing        (spacing between adject points, on both X and Y dimensions 1X2 vector)
%   .arenaSize      (optional). If not given the default is [96,32]
%  
% OUTPUT
% maskPos - NX2 matrix of positions in arena coordinates


if isfield(paramStruct, 'arenaSize')
    arenaSiz = paramStruct.arenaSize;
else
    arenaSiz  = [96, 32];
end

% indexing is off on purpose to convert spatial to matlab arena coordinates
xCrds = paramStruct.startPos(1):paramStruct.spacing(1):(paramStruct.startPos(1)+paramStruct.spacing(1)*paramStruct.gridSize(1)-1);
yCrds = paramStruct.startPos(2):paramStruct.spacing(2):(paramStruct.startPos(2)+paramStruct.spacing(2)*paramStruct.gridSize(2)-1);

if min(xCrds) < 1 || min(yCrds) < 1
    error('position outside arena range: lees than 1')
elseif max(xCrds) > arenaSiz(1)
    error('position outside arena range: exceed in X')
elseif max(yCrds) > arenaSiz(2)
    error('position outside arena range: exceed in Y')
end


[tempX, tempY] = meshgrid(xCrds, yCrds);

maskPos = [round(tempX(:)), round(tempY(:))];



end


