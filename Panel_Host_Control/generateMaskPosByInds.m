function structWmaskPosField  = generateMaskPosByInds(xRange, yRange)

% function structWmaskPosField  = generateMaskPosByInds(xrange, yrange)
%
% This function is designed to generate a structure with desired mask
% position base on the xrange, yrange inputs. It circumevents the usualy
% calculations for generating mask positions based on the stimulus itself
%
% INPUT
% 
% xRange(yRange) -          1XN(M) vector that includes the desired pixels indecies in the
%                           relveant dimension. 
%
% OUTPUT
% structWmaskPosField -     structure with maskPositions field that is
%                           generated based on specs. function outputs a structure so that it can be
%                           fed back into runPosFuncProtocol directly.
%
% NOTE! function has little input control since that should be dealt with
% in other functions

[xGrd, yGrd] = meshgrid(xRange, yRange);

structWmaskPosField.maskPositions = [xGrd(:), yGrd(:)];


end

