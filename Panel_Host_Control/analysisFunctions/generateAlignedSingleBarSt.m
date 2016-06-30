function singleBarSt = generateAlignedSingleBarSt(pStruct)

% function singleBarSt = generateAlignedSingleBarSt(pStruct)
%
% This function take the original protocolStruct from a SingleBarDiagCorr protocol 
% after appear and disappear have been added to the gratingTable and aligns
% it to bar appearance. Resulting structure is organized in position X stimDur structure
% Function also adds baseline subtracted data (baseline calculted by
% prestimBaseWin

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation
sdFac = 3;


% checking inputs
relTab = pStruct.gratingTable;

assert(ismember('appear', relTab.Properties.VariableNames), ...
           'gratingTable is missing appear variable')

assert(ismember('position', relTab.Properties.VariableNames), 'This function is designed for singleBar protocols only')

alignSt = alignProtocolDataByTable(pStruct, 'appear');

uPos = unique(relTab.position);
uDur = unique(relTab.stimDur);

singleBarSt = struct;
allBaseline = [];

for ii=1:length(uPos)
    
    for jj=1:length(uDur)
        
        relInd = ismember(relTab{:, {'position'; 'stimDur'}}, [uPos(ii), uDur(jj)], 'rows');
        
        singleBarSt(ii, jj).data = alignSt(relInd);
                
        relDat = singleBarSt(ii, jj).data.align.mean;
        relDatMed = singleBarSt(ii, jj).data.align.median;
        baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
        zeroInd = find(relDat(:,1) > 0, 1, 'first');
                
        baseVals = relDat(baseInds(1):baseInds(2), 2);
        baseSubResp = relDat;
        baseSubResp(:,2) = baseSubResp(:,2) - mean(baseVals);
        baseSubMed = relDatMed(:,2) - mean(baseVals);
                
        singleBarSt(ii, jj).subData.baseSub = baseSubResp;
        singleBarSt(ii, jj).subData.baseSubMed = baseSubMed;
        singleBarSt(ii, jj).subData.baseline = baseVals;
        singleBarSt(ii, jj).subData.zeroInd = zeroInd;
        singleBarSt(ii, jj).subData.length = size(baseSubResp, 1);
        
        allBaseline = vertcat(allBaseline, baseVals);
        
    end
    
end

baseSD = std(allBaseline);

% calculting max/min response

for ii=1:length(uPos)
    
    for jj=1:length(uDur)
        
        relDat = singleBarSt(ii,jj).subData;
        relInd = relDat.zeroInd;
        maxResp = quantile(relDat.baseSub(relInd:end, 2), 0.99);
        minResp = quantile(relDat.baseSub(relInd:end, 2), 0.01);
        
        if maxResp > baseSD * sdFac
            singleBarSt(ii,jj).maxResp = maxResp;
        else
            singleBarSt(ii,jj).maxResp = 0;
        end
        
        if abs(minResp) > baseSD * (sdFac-1) % since hyp is not sym to depol
            singleBarSt(ii,jj).minResp = minResp;
        else
            singleBarSt(ii,jj).minResp = 0;
        end
        
    end
    
end



end

