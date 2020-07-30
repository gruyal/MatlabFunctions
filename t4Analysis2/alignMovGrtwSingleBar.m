function alignMVPGrt = alignMovGrtwSingleBar(movGrtStruct, sbSt)

% function alignMVGrt = alignMovGrtwSingleBar(movGrtStruct, sbSt)
%
% this function is a modification of alignGratingPhasewSingleBar and alignMovGrtForAllProt
% combining them together since T4new does not have a gratingPhase
% protocol. 
%
% it also uses the new alignning functions (ending with 2) 
%
% INPUT
% movGrtStruct -        protocolStruct for moving grating
% sbST -                protocolStruct for singlebar 
%
% OUTPUT 

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

grtPC = movGrtStruct.inputParams.gridCenter; 
sbC = sbSt.inputParams.gridCenter; 
assert(unique(diff(vertcat(grtPC, sbC))) == 0, 'centers are not the same')

relO = movGrtStruct.inputParams.orientations; 

[alignSBSt, baseSB] = generateAlignedSingleBarStwMinMaxDiffWandV2(sbSt, 1);
relMaxExt = alignSBSt(end, end, end, end).maxExtPosVal;
relPD = sign(alignSBSt(end, end, end, end).minInhPosVal - relMaxExt);

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


relTab = movGrtStruct.gratingTable;

relTab.direction = relTab.direction * relPD; 
relTab.startPhase = relTab.startPhase * relPD; 

% adding norm positions to the table

grtPhase = (1:8)'; 
grtConvTab = table(grtPhase, normPosD, normPosB);

alignMVPGrt = struct; 
alignMVPGrt.gratingTable = relTab; 
alignMVPGrt.phaseTable = grtConvTab; 
alignMVPGrt.normPosWin = normPosWin; 

alignSt = alignProtocolDataByTable2(movGrtStruct, {'appear', 'appear', 'disappear'}, baseSB);


for ii=1:length(alignSt)

    alignMVPGrt.result(ii).data = alignSt(ii); 

    relDat = alignMVPGrt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = nanmean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    alignMVPGrt.result(ii).subData.baseSub = baseSubResp;
    alignMVPGrt.result(ii).subData.baseline = baseVal;
    alignMVPGrt.result(ii).subData.zeroInd = zeroInd;
    alignMVPGrt.result(ii).subData.length = size(baseSubResp, 1);

end



end
    

