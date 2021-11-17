function mvBarShiftSt = alignMovBarShiftwSingleBar(movBarShiftSt, sbSt)

% function combGrtStruct = alignMovBarShiftwSingleBar(movBarSt, grtMvSt, sbST)
%
% this functions takes the movimng bar data aligns it to the singlebar
% positions. 
% quick version from 2018 Oct 16
% will need to change when actually processing data
%
%
% INPUT
% movBarShiftSt -           moving bar shift protocol structure 
% sbSt -                    singlebar protocol structure
% 
% NOTE!         all protocols should have appear and disappear added to
%               their gratingTable
%
% OUTPUT
%
% TBD


pathFlag = 1; % 1 for T4 o for T5

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

mvPC = movBarShiftSt.inputParams.gridCenter; 
sbC = sbSt.inputParams.gridCenter;
sbO = sbSt.inputParams.orientations;

alignSBSt = generateAlignedSingleBarStwMinMaxDiffWandV(sbSt, pathFlag);
relMaxExt = alignSBSt(end, end, end, end).maxExtPosVal;
relPD = sign(alignSBSt(end, end, end, end).maxInhPosVal - relMaxExt);

% determining PD for movBar data

% if relPD > 0
%     pdO = sbO;
%     ndO = sbO+4;
% else
%     pdO = sbO+4;
%     ndO = sbO;
% end

assert(unique(diff(vertcat(mvPC, sbC))) == 0, 'centers are not the same')

% grating phase to singlebar converions (different for diagonal or
% non-diagonal orientation

mvsTab = movBarShiftSt.gratingTable; 
startPos = mvsTab.startPos; 
endPos = mvsTab.endPos; 
shiftPos = mvsTab.shiftPos; 
barWid = mvsTab.width;
uWid = unique(barWid); 
posSqLen = cellfun(@length, mvsTab.posSeq); 

for ww=1:length(uWid)
    tempSqLen = cellfun(@length, mvsTab.posSeq(mvsTab.width == uWid(ww))); 
    maxSeqPerWid(ww) = max(tempSqLen); 
end

if rem(sbO, 2) 
    diagF = 1;
else
    diagF = 0;
end

if relPD ==1
    normStPos = (startPos - relMaxExt) * relPD; 
    normEndPos = (endPos - relMaxExt) * relPD; 
    normShPos = (shiftPos - relMaxExt) * relPD; 
else
    normStPos = (startPos - relMaxExt - barWid+1) * relPD; 
    normEndPos = (endPos - relMaxExt - barWid+1) * relPD; 
    normShPos = (shiftPos - relMaxExt - barWid+1) * relPD; 
end
    


% adding diagonal flag

dirFlags = nan(length(startPos),2); 

for ii=1:length(startPos)
    
    relSqMax = maxSeqPerWid(uWid == barWid(ii)); 
    
    if posSqLen(ii) <= relSqMax/2 % only half trajectory
        dirFlags(ii, 1) = sign(normEndPos(ii) - normStPos(ii)); 
    else
        if isnan(normShPos(ii)) % when there is direct movement
            dirFlags(ii, :) = sign(normEndPos(ii) - normStPos(ii)); 
        else
            dirFlags(ii, 1) = sign(normShPos(ii) - normStPos(ii)); 
            dirFlags(ii, 2) = sign(normEndPos(ii) - normShPos(ii)); 
        end
    end
    
end


mvsTab = [mvsTab, table(normStPos, normEndPos, normShPos, dirFlags)];

alignSt = alignProtocolDataByTable(movBarShiftSt, 'appear');

mvBarShiftSt = struct; 

for ii=1:length(alignSt)
    
    mvBarShiftSt.result(ii).data = alignSt(ii); 
    
    tempT = movBarShiftSt.gratingTable(ii, :);
    nonShPI = ~cellfun(@(x) strcmp(x, 'shiftPos'), tempT.Properties.VariableNames); % removing shiftPos since it can be Nan and equal doesnt deal with it
    
    assert(isequal(alignSt(ii).table(1,nonShPI), tempT(1,nonShPI)), 'tables from alignSt and original struct are not equal')
    
    relDat = mvBarShiftSt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    mvBarShiftSt.result(ii).subData.baseSub = baseSubResp;
    mvBarShiftSt.result(ii).subData.baseline = baseVal;
    mvBarShiftSt.result(ii).subData.zeroInd = zeroInd;
    mvBarShiftSt.result(ii).subData.length = size(baseSubResp, 1);
    
end

mvBarShiftSt.table = mvsTab; 
mvBarShiftSt.diagFlag = diagF;




end

