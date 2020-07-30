function grtSt = alignGratingPhasewSingleBar(grtPhaseSt, sbSt)

% function combGrtStruct = alignGratingDatawSingleBar(grtPhaseSt, grtMvSt, sbST)
%
% this functions takes the grating phase data aligns it to the singlebar
% positions (based on which positions are ON or OFF and summs it to predict
% the moving grating responses 
%
%
% INPUT
% grtPhaseSt -              grating phase protocol structure 
% sbSt -                    singlebar protocol structure
% 
% NOTE!         all protocols should have appear and disappear added to
%               their gratingTable
%
% OUTPUT
%
% TBD

% verify center is the same


preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

grtPC = grtPhaseSt.inputParams.gridCenter; 
sbC = sbSt.inputParams.gridCenter; 

relO = grtPhaseSt.inputParams.orientations; 

alignSBSt = generateAlignedSingleBarStwMinMaxDiffWandV(sbSt);
relMaxExt = alignSBSt(end, end, end, end).maxExtPosVal;
relPD = sign(alignSBSt(end, end, end, end).maxInhPosVal - relMaxExt);

% % because of cell 15 - not worth it since not aligned with other
% if any(grtPC ~= sbC) 
%      offset = input('centers are not identical - add offset \n');
%      assert(round(offset) == offset, 'offset should be an integer')
% else
%      offset = 0;
% end

assert(unique(diff(vertcat(grtPC, sbC))) == 0, 'centers are not the same')

% grating phase to singlebar converions (different for diagonal or
% non-diagonal orientation

grtPhase = (1:8)';
if rem(relO, 2) 
    
    sbPosD = {[-6:-4, 1:4]; [-6:-3, 2:5]; [-5:-2, 3:6]; [-4:-1, 4:6]; [-3:0, 5:6]; [-2:1, 6]; [-1:2, -6]; [0:3, -6:-5]}; 
    sbPosB = cellfun(@(x) setdiff(-6:6, x), sbPosD, 'uniformoutput', 0); 
    posWin = -6:6;
else
    
    sbPosD = {[-4, 1:4]; [-4:-3, 2:4]; [-4:-2, 3:4]; [-4:-1, 4]; -3:0; -2:1; -1:2; 0:3}; % dark positions
    sbPosB = cellfun(@(x) setdiff(-4:4, x), sbPosD, 'uniformoutput', 0); % bright positions    
    posWin = -4:4;
end


normPosD = cellfun(@(x) (x - relMaxExt) * relPD, sbPosD, 'uniformoutput', 0);
normPosB = cellfun(@(x) (x - relMaxExt) * relPD, sbPosB, 'uniformoutput', 0);

normPosWin = (posWin - relMaxExt) * relPD; 

% adding grtComb if it is absent
grtTab = grtPhaseSt.gratingTable;

if ~ismember('grtComb', grtTab.Properties.VariableNames)
    tempComb = round(2*(grtTab.FBval + grtTab.SBVal));
    tempComb(tempComb == 2) = 4; 
    tempComb(tempComb == 3) = 2; 
    tempComb(tempComb == 4) = 3;
    
    grtComb = tempComb; 
    grtTab = [grtTab, table(grtComb)];
    
end
    

% adding norm positions to the table
grtConvTab = table(normPosD, normPosB);
grtTab = [grtTab, grtConvTab(grtTab.phase, :)];

alignSt = alignProtocolDataByTable(grtPhaseSt, 'appear');

grtSt = struct; 

for ii=1:length(alignSt)
    
    grtSt.result(ii).data = alignSt(ii); 
    
    tempT = grtPhaseSt.gratingTable(ii, :);
    
    assert(isequal(alignSt(ii).table, tempT), 'tables from alignSt and original struct are not equal')
    
    relDat = grtSt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    grtSt.result(ii).subData.baseSub = baseSubResp;
    grtSt.result(ii).subData.baseline = baseVal;
    grtSt.result(ii).subData.zeroInd = zeroInd;
    grtSt.result(ii).subData.length = size(baseSubResp, 1);
    
end

grtSt.table = grtTab; 
grtSt.normPosWin = normPosWin; 


end

