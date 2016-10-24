function paddeVec = padRespVecGen(respSt, newZeroInd, newTotLength)

% function paddeVec = padRespVecGen(respSt, newZeroInd, newTotLength)
%
% This is a general modification of padRespVec. It requires a simpler
% structure as input but functions exactly the same. 
% 
%
% it takes a response structure and moves its zeroInd to the newZeroInd
% and makes its total length newTotLength.
%
% INPUT
% 
% respSt -          Should have the following fields
%   .data -         NX1 vector to shift to the required index
%   .zeroInd -      integer. index the spcifies the old zeroIndex
% 
% newZeroInd -      integer. the index to which the old zero will be
%                   aligned to (either by prepadding zeros or clipping)
% newTotLen -       integer. new total length 
%
% OUTPUT
% 
% paddedVec -       just the resp itself < respSt.subData.baseSub(:,2)>
%                   after it had been modified/shifted



oldZeroInd = respSt.zeroInd;
relSt = respSt.data;



preDiff = newZeroInd - oldZeroInd;
                
if preDiff > 0
    paddeVec = padarray(relSt, [preDiff, 0], 0, 'pre');
else
    paddeVec = relSt(abs(preDiff)+1:end);
end

oldLen = length(relSt) + preDiff;

postDiff = newTotLength - oldLen;

if postDiff > 0
    paddeVec = padarray(paddeVec, [postDiff, 0], 0, 'post');
else
    paddeVec = paddeVec(1:end-abs(postDiff));
end




end