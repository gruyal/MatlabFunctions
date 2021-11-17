function mvBarSt = alignMovBarwSingleBar(movBarSt, sbSt)

% function combGrtStruct = alignGratingDatawSingleBar(movBarSt, grtMvSt, sbST)
%
% this functions takes the movimng bar data aligns it to the singlebar
% positions. 
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

alignSBSt = generateAlignedSingleBarStwMinMaxDiffWandV(sbSt);
relMaxExt = alignSBSt(end, end, end, end).maxExtPosVal;
relPD = sign(alignSBSt(end, end, end, end).minInhPosVal - relMaxExt);

% determining PD for movBar data

if relPD > 0
    pdO = sbO;
    ndO = sbO+4;
else
    pdO = sbO+4;
    ndO = sbO;
end

assert(unique(diff(vertcat(mvPC, sbC))) == 0, 'centers are not the same')

% grating phase to singlebar converions (different for diagonal or
% non-diagonal orientation

if rem(sbO, 2) 
    diagF = 1;
    posWin = -6:6;
else
    diagF = 0;
    posWin = -4:4;
end

normPosWin = (posWin - relMaxExt) * relPD; 

% adding diagonal flag
mvTab = movBarSt.gratingTable;
PDFlag = mvTab.orient == pdO;
NDFlag = mvTab.orient == ndO;

mvTab = [mvTab, table(PDFlag, NDFlag)];

alignSt = alignProtocolDataByTable(movBarSt, 'appear');

mvBarSt = struct; 

for ii=1:length(alignSt)
    
    mvBarSt.result(ii).data = alignSt(ii); 
    
    tempT = movBarSt.gratingTable(ii, :);
    
    assert(isequal(alignSt(ii).table, tempT), 'tables from alignSt and original struct are not equal')
    
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
    
end

mvBarSt.table = mvTab; 
mvBarSt.normPosWin = normPosWin;
mvBarSt.posWin  = posWin; 
mvBarSt.relMaxExt = relMaxExt;
mvBarSt.relPD = relPD;
mvBarSt.centerNormPos = normPosWin(ceil(length(normPosWin)/2));
mvBarSt.diagFlag = diagF;




end

