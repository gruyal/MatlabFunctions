function corrPStruct = correctCorruptTimeStim(pStruct, stimNum, refRep)

% function correctCorruptTimeStim(pStruct, stimNum, refStim)
%
% This function deals is meant to correct for stim in which it is clear the
% response is shifted in time due to a corrupt time vector. It should be
% used before the alignProtocolDataByTable2 (though alignning to find stim
% to correct is useful - since can plot all with plotAlignedResponsesWRepeats)
%
% INPUT
%
% pStruct -         protocolStruct with the stim to correct
% stimNum -         number for the stimulus to correct (from gratingTable)
% refRep -          stimulus repeat to use as reference (againt which cross
%                   correlation would be calculated)
%                   (should be given as index for repeats number and not
%                   protocolStruct stim number)
%                   If relRep is given as negative number that repeat would
%                   be removed - added for grating stim that dont make
%                   sense (look like another stim)
%
% OUTPUT
% corrPStruct -     protocolStruct with the corrected (shifted repeats) in
%                   the proper positions

corrPStruct = pStruct; 

stimIndSt = getStimInds(pStruct, [stimNum, nan, nan, nan]);
relInds = stimIndSt.inds; 

assert(abs(refRep) <= length(relInds), 'refRep is bigger than number of repeats to this stimulus')


if refRep > 0 % aligns to that relRep 
    baseD = pStruct.stim(relInds(refRep)).data{1}(:,2);
    meanBase = mean(baseD); 
    sdBase = std(baseD); 
    zBase = (baseD - meanBase) ./ sdBase; 

    for ii=1:length(relInds)

        tempD = pStruct.stim(relInds(ii)).data{1}(:,2);
        meanTemp = mean(tempD); 
        sdTemp = std(tempD); 
        zTemp = (tempD - meanTemp) ./ sdTemp; 

        [xCorrRes, xCorrLag] = xcorr(zBase, zTemp); 
        [~, corrMaxI] = max(xCorrRes); 
        diffSamp = xCorrLag(corrMaxI);

        switch sign(diffSamp)

            case 0 % no change
                corrPStruct.stim(relInds(ii)).data{1}(:,2) = tempD; % no change 

            case -1
                corrPStruct.stim(relInds(ii)).data{1}(:,2) = [tempD(abs(diffSamp):end); nan(abs(diffSamp)-1,1)]; % move forward

            case 1
                corrPStruct.stim(relInds(ii)).data{1}(:,2) = [nan(diffSamp,1); tempD(1:end-diffSamp)]; % move backward
        end

    end
elseif refRep < 0 % removes the relRep (by changing the baseline so that it will be weeded out 
    
    tempD = pStruct.stim(relInds(abs(refRep))).data{1}(:,2);
    apPos = pStruct.gratingTable.appear(stimNum);
    apPosInd = pStruct.stim(relInds(abs(refRep))).data{2}(:,2) == apPos;
    timeInd = pStruct.stim(relInds(abs(refRep))).data{2}(apPosInd,1);
    datInd = pStruct.stim(relInds(abs(refRep))).data{1}(:,1) < timeInd; 
    tempD(datInd) = 0; % will exclude this repeat in alignProtocolDataByTable2
    corrPStruct.stim(relInds(abs(refRep))).data{1}(:,2) = tempD; 
    
end
    







end