function alignMVGrt = alignMovGrtForAllProt(movGrtStruct)

% function alignMVGrt = alignMovGrtForAllProt(movGrtStruct)
%
% this function is used to aligned a moving grating protocol for the
% allProtocol structure. it is not meant for comparison between cells but
% within a cell (since the grating will be aligned to the appear frame)

relTab = movGrtStruct.gratingTable;
    
if ~ismember('grtComb', relTab.Properties.VariableNames)
    tempComb = round(2*(relTab.FBval + relTab.SBVal));
    tempComb(tempComb == 2) = 4; 
    tempComb(tempComb == 3) = 2; 
    tempComb(tempComb == 4) = 3;

    grtComb = tempComb; 
    relTab = [relTab, table(grtComb)];
end

alignMVGrt = struct; 
alignMVGrt.gratingTable = relTab; 

alignSt = alignProtocolDataByTable(movGrtStruct, 'appear');

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

for ii=1:length(alignSt)

    grtSt.result(ii).data = alignSt(ii); 

    relDat = grtSt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    alignMVGrt.result(ii).subData.baseSub = baseSubResp;
    alignMVGrt.result(ii).subData.baseline = baseVal;
    alignMVGrt.result(ii).subData.zeroInd = zeroInd;
    alignMVGrt.result(ii).subData.length = size(baseSubResp, 1);

end



end