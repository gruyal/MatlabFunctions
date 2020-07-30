function meanMaxMADTab = estimateMeanvsSSEperStim(allProtStruct)

% function meanVsSseAll = estimateMeanvsSSEperStim(allProtStruct)
% 
% this function takes the allProtSt (created using
% combiningAllProtStructScript) and calculates the mean square error
% between the mean of each response and the valid (non excluded) repeats
% that constitute the mean. The calculation is only performed in the max
% response (since it is to be used the max plot for modelling comparison) 
%
% NOTE !!!      window is calculated for T5 - for T4 flip valFactors (longer for
%               the dark stim) 
%
% INPUT
%
% allProtStruct -       structure containing all the protocol for a
%                       particulat cell (created w/ combiningAllProtStructScript). 
% OUTPUT
% 
% meanMaxMADTab -       Table containg stimulus type, index for each
%                       stimulus type (linking back to original stimulus
%                       table), meanMax resp, baseline, baseline subtracted meanMax, 
%                       repeats max (see note 2), and mean absolute
%                       deviation. 
%
% NOTE 2 -              repeats max are calculated around the time point in
%                       which max response is found for the mean (not using the same window)
%                       this is done to facilitate a better comparison
%                       between mean and repeats (since they are already
%                       aligned)


% protInd 
%   SB  1
%   MB  2
%   MM  3
%   SG  4


% calculating error in sb

numSBStim = height(allProtStruct.sbTable); 
Index = allProtStruct.sbTable.index; 
sbTab = allProtStruct.sbTable; 
unitStep = 20; %in ms

baseMMVals = 0:5:25; 

winSt = 40; % in ms 

% multipliers for max windows


sbValFac = [2,4]; % multiplier for window (for T5)
mbValFac = [1, 1.2];
mmFac = [1.5, 2];
repWinBuff = 50; % in ms
relGrtComb = 2; % change to 1 for T4

MeanMax = zeros(numSBStim, 1); 
BaseSubMeanMax = zeros(numSBStim, 1); 
Baseline = zeros(numSBStim, 1); 
RepMax = cell(numSBStim, 1); 
MAD = zeros(numSBStim, 1); 

for ii=1:numSBStim
    
    % since for SB not all the structure is filled
    sbInds = sbTab.inds(ii,:);
    relResSt = allProtStruct.sbResult(sbInds(1), sbInds(2), sbInds(3), sbInds(4)); 
    relDat = relResSt.data.align;
    exReps = relResSt.data.exclude; 
    relBase = mean(relResSt.subData.baseline); 
    Baseline(ii) = relBase; 
    %used to detemine window for max
    relDur = sbTab.stimDur(ii); 
    valInd = sbTab.value(ii) + 1; % since value is 0 and 1
    
    relWin  = [winSt, relDur*1000*sbValFac(valInd)+winSt];

    relMean = relDat.mean; 
    relInds = arrayfun(@(x) find(relMean(:,1) > x, 1, 'first'), relWin);
    [meanMax, mmInd] = max(relDat.mean(relInds(1):relInds(2), 2)); 
    peakTime = relDat.mean(relInds(1) + mmInd, 1);
    MeanMax(ii) = meanMax; 
    BaseSubMeanMax(ii) = meanMax - relBase; 
    
    allReps = 1:length(relDat.rep);
    relReps = setdiff(allReps, exReps); 
    numReps = length(relReps);
    tempMaxVec = zeros(1, numReps); 
    
    for jj=1:numReps
        
        tempRep = relDat.rep(relReps(jj)).data; 
        repInds = arrayfun(@(x) find(tempRep(:,1) > x, 1, 'first'), [peakTime - repWinBuff, peakTime + repWinBuff]);
        
        tempMax = max(tempRep(repInds(1):repInds(2), 2)); 
        
        tempMaxVec(jj) = tempMax; 
        
    end
    
    RepMax{ii} = tempMaxVec; 
    MAD(ii) = sum(abs(tempMaxVec - meanMax)) / numReps; 
    
end
        
ProtType = repmat({'SB'}, numSBStim, 1); 
ProtInd = ones(numSBStim, 1) * 1; 
sbMADTab = table(ProtType, ProtInd, Index, Baseline, MeanMax, BaseSubMeanMax, RepMax, MAD); 


% moving bar calculation


numMBStim = height(allProtStruct.mvTable); 
Index = allProtStruct.mvTable.index;

MeanMax = zeros(numMBStim, 1); 
BaseSubMeanMax = zeros(numMBStim, 1); 
Baseline = zeros(numMBStim, 1); 
RepMax = cell(numMBStim, 1); 
MAD = zeros(numMBStim, 1); 

for ii=1:numMBStim
    
    relResSt = allProtStruct.mvResult(ii); 
    
    relDat = relResSt.data.align;
    exReps = relResSt.data.exclude; 
    relBase = relResSt.subData.baseline; 
    Baseline(ii) = relBase; 
    
    relTab = relResSt.data.table; 
    stimTotDur = ((relTab.disappear - relTab.appear) / relTab.framePerStep) * relTab.stepDur * 1000;
    valInd = relTab.value +1; 
    stimEnd = stimTotDur * mbValFac(valInd) + winSt; 
    
    relWin  = [winSt, stimEnd];
    
    relMean = relDat.mean; 
    relInds = arrayfun(@(x) find(relMean(:,1) > x, 1, 'first'), relWin, 'uniformoutput', 0);
    % in case the window is byond the end
    if isempty(relInds{2})
        relInds = [relInds{1}, size(relMean, 1) - 500]; % almost towards the end
    else
        relInds = [relInds{1}, relInds{2}];
    end
    [meanMax, mmInd] = max(relDat.mean(relInds(1):relInds(2), 2)); 
    peakTime = relDat.mean(relInds(1) + mmInd);  
    MeanMax(ii) = meanMax; 
    BaseSubMeanMax(ii) = meanMax - relBase; 
    
    allReps = 1:length(relDat.rep);
    relReps = setdiff(allReps, exReps); 
    numReps = length(relReps);
    tempMaxVec = zeros(1, numReps); 
    
    for jj=1:numReps
        
        tempRep = relDat.rep(relReps(jj)).data; 
        if peakTime + repWinBuff > tempRep(end,1)
            repWin  = [peakTime - repWinBuff, tempRep(end,1) - 20]; % in ms
        else
            repWin = [peakTime - repWinBuff, peakTime + repWinBuff];
        end
        
        repInds = arrayfun(@(x) find(tempRep(:,1) > x, 1, 'first'), repWin);
        
        tempMax = max(tempRep(repInds(1):repInds(2), 2)); 
        
        tempMaxVec(jj) = tempMax; 
        
    end
    
    RepMax{ii} = tempMaxVec; 
    MAD(ii) = sum(abs(tempMaxVec - meanMax)) / numReps; 
    
end
    
ProtType = repmat({'MB'}, numMBStim, 1);
ProtInd = ones(numMBStim, 1) *2; 
mbMADTab = table(ProtType, ProtInd, Index, Baseline, MeanMax, BaseSubMeanMax, RepMax, MAD); 


% minMot calculation


numMMStim = height(allProtStruct.mmTable); 
Index = allProtStruct.mmTable.index;

MeanMax = zeros(numMMStim, 1); 
BaseSubMeanMax = zeros(numMMStim, 1); 
Baseline = zeros(numMMStim, 1); 
RepMax = cell(numMMStim, 1); 
MAD = zeros(numMMStim, 1); 

for ii=1:numMMStim
    
    relResSt = allProtStruct.mmResult(ii); 
    
    relDat = relResSt.data.align;
    exReps = relResSt.data.exclude; 
    relBase = relResSt.subData.baseline; 
    Baseline(ii) = relBase; 
    
    relTab = relResSt.data.table; 
    valInd = relTab.FBVal + 1; 
    stimEnd = (relTab.sDisappear - relTab.fAppear) * unitStep * mmFac(valInd) + winSt;
    
    relWin  = [winSt, stimEnd];
    
    relMean = relDat.mean; 
    relInds = arrayfun(@(x) find(relMean(:,1) > x, 1, 'first'), relWin, 'uniformoutput', 0);
    % in case the window is byond the end
    if isempty(relInds{2})
        relInds = [relInds{1}, size(relMean, 1) - 500]; % almost towards the end (in samples)
    else
        relInds = [relInds{1}, relInds{2}];
    end
    [meanMax, mmInd] = max(relDat.mean(relInds(1):relInds(2), 2)); 
    peakTime = relDat.mean(relInds(1) + mmInd);  
    MeanMax(ii) = meanMax; 
    BaseSubMeanMax(ii) = meanMax - relBase; 
    
    allReps = 1:length(relDat.rep);
    relReps = setdiff(allReps, exReps); 
    numReps = length(relReps);
    tempMaxVec = zeros(1, numReps); 
    
    for jj=1:numReps
        
        tempRep = relDat.rep(relReps(jj)).data; 
        if peakTime + repWinBuff > tempRep(end,1)
            repWin  = [peakTime - repWinBuff, tempRep(end,1) - 20]; % in ms
        else
            repWin = [peakTime - repWinBuff, peakTime + repWinBuff];
        end
        repInds = arrayfun(@(x) find(tempRep(:,1) > x, 1, 'first'), repWin);
        
        tempMax = max(tempRep(repInds(1):repInds(2), 2)); 
        
        tempMaxVec(jj) = tempMax; 
        
    end
    
    RepMax{ii} = tempMaxVec; 
    MAD(ii) = sum(abs(tempMaxVec - meanMax)) / numReps; 
    
end
    
ProtType = repmat({'MM'}, numMMStim, 1); 
ProtInd = ones(numMMStim, 1) *3; 
mmMADTab = table(ProtType, ProtInd, Index, Baseline, MeanMax, BaseSubMeanMax, RepMax, MAD); 



% grating phase calculation


if isfield(allProtStruct, 'gpResult')

    numGPStim = height(allProtStruct.gpTable); 
    gpTab = allProtStruct.gpTable; 
    Index = allProtStruct.gpTable.index;

    MeanMax = zeros(numGPStim, 1); 
    BaseSubMeanMax = zeros(numGPStim, 1); 
    Baseline = zeros(numGPStim, 1); 
    RepMax = cell(numGPStim, 1); 
    MAD = zeros(numGPStim, 1); 

    for ii=1:numGPStim

        relResSt = allProtStruct.gpResult(ii); 

        relDat = relResSt.data.align;
        exReps = relResSt.data.exclude; 
        relBase = relResSt.subData.baseline; 
        Baseline(ii) = relBase; 

        relTab = gpTab(ii,:);
        stimTotDur = (relTab.disappear - relTab.appear) * unitStep;

        if relTab.grtComb == relGrtComb 
            stimEnd = stimTotDur + 200; % in ms
        else
            stimEnd = stimTotDur + 20;
        end

        relWin  = [winSt, stimEnd];

        relMean = relDat.mean; 
        relInds = arrayfun(@(x) find(relMean(:,1) > x, 1, 'first'), relWin);
        [meanMax, mmInd] = max(relDat.mean(relInds(1):relInds(2), 2)); 
        peakTime = relDat.mean(relInds(1) + mmInd);  
        MeanMax(ii) = meanMax; 
        BaseSubMeanMax(ii) = meanMax - relBase; 

        allReps = 1:length(relDat.rep);
        relReps = setdiff(allReps, exReps); 
        numReps = length(relReps);
        tempMaxVec = zeros(1, numReps); 

        for jj=1:numReps

            tempRep = relDat.rep(relReps(jj)).data; 
            repInds = arrayfun(@(x) find(tempRep(:,1) > x, 1, 'first'), [peakTime - repWinBuff, peakTime + repWinBuff]);

            tempMax = max(tempRep(repInds(1):repInds(2), 2)); 

            tempMaxVec(jj) = tempMax; 

        end

        RepMax{ii} = tempMaxVec; 
        MAD(ii) = sum(abs(tempMaxVec - meanMax)) / numReps; 

    end

    ProtType = repmat({'GP'}, numGPStim, 1); 
    ProtInd = ones(numGPStim, 1) *4; 
    gpMADTab = table(ProtType, ProtInd, Index, Baseline, MeanMax, BaseSubMeanMax, RepMax, MAD); 
else
    gpMADTab = [];
end

meanMaxMADTab = [sbMADTab; mbMADTab; mmMADTab; gpMADTab];

grpByMM = nan(height(meanMaxMADTab), 1); 

for ii=1:length(baseMMVals)-1
    tmpBI = meanMaxMADTab.BaseSubMeanMax > baseMMVals(ii) & meanMaxMADTab.BaseSubMeanMax <= baseMMVals(ii+1);
    grpByMM(tmpBI) = ii;
end

meanMaxMADTab = [meanMaxMADTab, table(grpByMM)];


end
