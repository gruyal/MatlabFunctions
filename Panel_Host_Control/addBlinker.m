function modFramesMat = addBlinker(framesMat, gsLevel)

% function modFramesMat = addBlinker(framesMat, gsLevel)
%
% This function adds a blinker to the framesMat (toggles each frame)
%
% INPUT
%
% framesMat -       32X96XN matrix (frames to be presented)
% gsLevel -         gray Scale level for framesMat.
%
% OUTPUT
%
% modFramesMat -    modified framesMat, with MIN/MAX pixels added @ x=1:2,
%                   y=20:22 at odd/even frames respectively


modFramesMat = framesMat;
maxVal = 2^gsLevel-1;


modFramesMat(1:2, 20:22, 1:2:end) = 0;
modFramesMat(1:2, 20:22, 2:2:end) = maxVal;







end