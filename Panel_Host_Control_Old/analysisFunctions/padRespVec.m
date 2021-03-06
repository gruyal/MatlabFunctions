function paddedVec = padRespVec(respSt, newZeroInd, newTotLength)

% function paddeVec = padRespVec(respSt, newZeroInd, newTotLength)
%
% This function is to be used internally to generate linear comparison
% between bar flashes and movingbar/minMot stim
%
% it takes a response structure and moves its zeroInd to the newZeroInd
% and makes its total length newTotLength.
%
% INPUT
% 
% respSt -          generated by passing a protocol (singleBar) through generateAlignedSingleBarSt
%                   but choosing just one position and duration
% newZeroInd -      integer. the index to which the old zero will be
%                   aligned to (either by prepadding zeros or clipping)
% newTotLen -       integer. new total length 
%
% OUTPUT
% 
% paddedVec -       just the resp itself < respSt.subData.baseSub(:,2)>
%                   after it had been modified/shifted


% changed the location of zeroInd in some function, so included a
% modification to look for it at 2 separate layers

if isfield(respSt, 'subData')
    oldZeroInd = respSt.subData.zeroInd; 
    relSt = respSt.subData;
elseif isfield(respSt, 'zeroInd')
    oldZeroInd = respSt.zeroInd;
    relSt = respSt;
else
    error('could not locate original zeroInd')
end


preDiff = newZeroInd - oldZeroInd;
                
if preDiff > 0
    paddedVec = padarray(relSt.baseSub(:,2), [preDiff, 0], 0, 'pre');
else
    paddedVec = relSt.baseSub(abs(preDiff)+1:end,2);
end

oldLen = relSt.length + preDiff;

postDiff = newTotLength - oldLen;

if postDiff > 0
    paddedVec = padarray(paddedVec, [postDiff, 0], 0, 'post');
else
    paddedVec = paddedVec(1:end-abs(postDiff));
end

paddedVec = extrapolateShiftedVec(paddedVec);


end