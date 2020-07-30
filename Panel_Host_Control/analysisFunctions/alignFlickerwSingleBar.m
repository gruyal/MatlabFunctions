function flkBarSt = alignFlickerwSingleBar(flickerSt, sbSt)

% function combGrtStruct = alignGratingDatawSingleBar(movBarSt, grtMvSt, sbST)
%
% this functions takes the movimng bart data aligns it to the singlebar
% positions. 
%
%
% INPUT
% flickerSt -               flicker protocol structure 
% sbSt -                    singlebar protocol structure
% 
% NOTE!         all protocols should have appear and disappear added to
%               their gratingTable
%
% OUTPUT
%
% TBD



preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

flkPC = flickerSt.inputParams.gridCenter; 
sbC = sbSt.inputParams.gridCenter;
sbO = sbSt.inputParams.orientations;

assert(unique(diff(vertcat(flkPC, sbC))) == 0, 'centers are not the same')

if rem(sbO, 2) 
    diagF = 1;
else
    diagF = 0;
end


alignSBSt = generateAlignedSingleBarStwMinMaxDiffWandV(sbSt);
relMaxExt = alignSBSt(end, end, end, end).maxExtPosVal;
relPD = sign(alignSBSt(end, end, end, end).minInhPosVal - relMaxExt);

% adding corrected positions and the first frame of OFF appearing 
flkTab = flickerSt.gratingTable;
origPos = flkTab.position; 
barWid = flkTab.width; 

normPos = (origPos - relMaxExt - barWid + 1) * relPD; 

% since cell 1 starts with an off bar (rest with on)
offApp = flkTab.offAppear; 
firstOff = cellfun(@(x) x(1), offApp); 

flkTab = [flkTab, table(normPos, firstOff)];

flickerSt.gratingTable = flkTab; 

alignSt = alignProtocolDataByTable(flickerSt, 'firstOff');

flkBarSt = struct; 

for ii=1:length(alignSt)
    
    flkBarSt.result(ii).data = alignSt(ii); 
    
    tempT = flickerSt.gratingTable(ii, :);
    
    assert(isequal(alignSt(ii).table, tempT), 'tables from alignSt and original struct are not equal')
    
    relDat = flkBarSt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    flkBarSt.result(ii).subData.baseSub = baseSubResp;
    flkBarSt.result(ii).subData.baseline = baseVal;
    flkBarSt.result(ii).subData.zeroInd = zeroInd;
    flkBarSt.result(ii).subData.length = size(baseSubResp, 1);
    
end

flkBarSt.table = flkTab; 
flkBarSt.uNormPos = unique(normPos);
flkBarSt.diagFlag = diagF;




end

