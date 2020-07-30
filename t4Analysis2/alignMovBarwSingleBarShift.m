function mvBarSt = alignMovBarwSingleBarShift(movBarSt, sbSt)

% function mvBarSt = alignMovBarwSingleBarShift(movBarSt, sbST)
%
% this functions takes the moving bar shift data aligns it to the singlebar
% positions. (it is a modification of alignMovBarwSingleBar
% 
%  !!! NOTE currently only dealing w/ full trajectories and not with their
%  components
%
%
% INPUT
% movBarSt -                moving bar protocol structure 
% sbSt -                    singlebar protocol structure
% 
% NOTE!         all protocols should have appear and disappear added to
%               their gratingTable
%
% OUTPUT
%
% TBD



preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

mvPC = movBarSt.inputParams.gridCenter; 
sbC = sbSt.inputParams.gridCenter;
sbO = sbSt.inputParams.orientations;

[alignSBSt, baseSt] = generateAlignedSingleBarStwMinMaxDiffWandV2(sbSt, 1);
relMaxExt = alignSBSt(end, end, end, end).maxExtPosVal;
relPD = sign(alignSBSt(end, end, end, end).minInhPosVal - relMaxExt);

if rem(sbO, 2) 
    diagF = 1;
else
    diagF = 0;
end


assert(unique(diff(vertcat(mvPC, sbC))) == 0, 'centers are not the same')

% grating phase to singlebar converions (different for diagonal or
% non-diagonal orientation

movBarSt.gratingTable.shiftPos(isnan(movBarSt.gratingTable.shiftPos)) = 999; 

alignSt = alignProtocolDataByTable2(movBarSt, {'appear', 'appear', 'disappear'}, baseSt);

mvTab = movBarSt.gratingTable;
newMVSTab = normMovBarShiftTable(mvTab, relPD, relMaxExt); 
mvBarSt = struct; 

for ii=1:length(alignSt)
    
    mvBarSt.result(ii).data = alignSt(ii); 
    
    tempT = mvTab(ii, :);
    
    assert(isequal(alignSt(ii).table, tempT), 'tables from alignSt and original struct are not equal in row %d', ii)
    
    relDat = mvBarSt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    mvBarSt.result(ii).subData.baseSub = baseSubResp;
    mvBarSt.result(ii).subData.baseline = baseVal;
    mvBarSt.result(ii).subData.zeroInd = zeroInd;
    mvBarSt.result(ii).subData.length = size(baseSubResp, 1);
    mvBarSt.result(ii).table = tempT; % just as sanity check - since in some cases I am excluding stim 
end

mvBarSt.table = newMVSTab; 
mvBarSt.oldTable = mvTab; 
mvBarSt.diagFlag = diagF;




end

