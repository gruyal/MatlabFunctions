function [alignStruct, alignTable, varargout] = alignProtocolDataByTableG4(pStruct, relVarNames, baseSDStruct)

% function alignStruct = alignProtocolDataByTableG4(pStruct, relVarName, baseSDStruct)
%
% This function uses getAlignedStimDataByTable2 to align the entire
% protocolStruct and outputs an aligned stucutre. it was modified from v2 
% to include refrence to maskPos and maskOrt also (needed for H2
% recordings)
%
% This is a modified version of the function designed to work with
% getAlignedStimDataByTable2 which reports also the postStim baseline. 
% This function also uses a different calcluation of which repeats to
% exclude base on both pre and post baselines
%
%
% INPUT
%
% pStruct -             protocolStruct w/.stim, .data and .gratingTable in it
% relVarName -          variable name in according to which data will be
%                       aligned. In version 2 has to be 3 vars long 
%                       1 - posVal to align to 
%                       2 - posVal to calculate preBaseline from 
%                       3 - posVal to calculate postBaseline from (after
%                       timeBuff delay) (see getAlignedStimDataByTable)
%       Note!!      Name should come from gratingTable variable and that variable should have 
%                   a relevant controller frame number (starting from 0) accrding to which data would be aligned
% baseSDStruct -        (optional) baseThresh struct from running this function on a
%                       different protocol from the same cell (in case this
%                       cell doesnt have enough stimuli to properly
%                       estimate the baseline SD). The final SD would a
%                       weighted sum of both SDs with the weights derived
%                       from the adjusted R square of the fits. 
%       NOTE: designed to be used with singleBar that has most stim 
%
%       NOTE: baseSDStruct would be used only if difference between low and
%       high is higher than 35mV (which means SD estimate is poor) 
%
%
% Removed:
% meanTreshFlag -       (optional). logical. If TRUE uses the high value
%                       for response treshold. if not lower. The different values are aim the
%                       control for stimuli that are long (moving bar) versus short pulsed type
%                       of stim (flicker or singleBar). Default 0. 
%
% OUTPUT
%
% alignStruct -         structure with the same fields as output of
%                       getAlignedDataByTable plus the relevant table row
%
% Note added medFindSpkPar to alignStruct as these are the filtering
%      parameters for medFiltAndFindSpikes
%
% alignTab -            Added since in G4 maskPos and maskOrt are sometimes
%                       used. This table integrates gratingTable with
%                       maskOrt and maskPos (relInds 3 and 4) to create a
%                       golbal table that is a reference into alignStruct
%
% stimToCorrInds -      indices for stim that were deemed inadmissable


if nargin < 3
    for ii=1:2
        baseSDStruct(ii).low = -65;
        baseSDStruct(ii).high = -65;
        baseSDStruct(ii).adjRSq = 0; % weight for this would be zero
    end
end


parSt.medWinDef = 351;
parSt.threshVDef = 4.5;
parSt.smWin = 21;

sigFac = 3.5; % factor to multifly sigma and exclude everything above mu + sigFac * sigma (from pre and post baseline calculation)
sigRng = 4.5; % used in the same claculation if the var is low (will allow baseline mean +/- this)

assert(isfield(pStruct, 'gratingTable'), 'pStruct is missing gratingTable field')

stimTable = pStruct.gratingTable;
stimOrt = pStruct.orientations;
stimMaskPos = pStruct.maskPositions;

assert(all(ismember(relVarNames, stimTable.Properties.VariableNames)), 'gratingTable does not contain one or all of relVarNames')
assert(length(relVarNames) == 3, 'relVarNames must be of length 3 - align, start and end names (align and start can be the same')
assert(ismember('index', stimTable.Properties.VariableNames), 'gratingTable does not contain index variable')

alignStruct = struct;
alignTable = table;

newIndex = 0;
for ii=1:height(stimTable)
    
    relTable = stimTable(ii, :);
    relVal = relTable{:, relVarNames};
    
    for oo=1:length(stimOrt)
        
        relOrt = stimOrt(oo);
        
        for pp=1:size(stimMaskPos,1)
            
            newIndex = newIndex+1;
            relPos = stimMaskPos(pp,:);
            relInd = [relTable.index, nan, oo, pp];
            
            alignStruct(newIndex).align = getAlignedStimDataByTableG4(pStruct, relInd, relVal, parSt);
            relRow = [table(newIndex), relTable, table(relOrt, relPos, relInd)];
            alignStruct(newIndex).table = relRow; 
            alignStruct(newIndex).exclude = [];
            alignStruct(newIndex).numReps = length(alignStruct(newIndex).align.rep);
            alignStruct(newIndex).medFindSpkPar = parSt;
            alignTable = [alignTable; relRow];
        end
        
    end
    
end

totReps = sum([alignStruct.numReps]); 

numEdges = floor(totReps/20 / 10) * 10; % worked for SB that has a lot of repeats 

if numEdges < 12 % since above calculation is for protocols with totReps > 500
    numEdges = floor(totReps/10); 
    if numEdges < 5
        numEdges = 5; % so gaussian could be fit
        warning('numEdges too low, increased to 5')
    end
        
end
    

% excluding repeats that were too noisy

allPreMean = nan(1, totReps);
allPostMean = nan(1, totReps);

count = 0;
for ii=1:length(alignStruct)
    tempRepNum = length(alignStruct(ii).align.rep); 
    for jj=1:tempRepNum
        count = count+1;
        allPreMean(count) = alignStruct(ii).align.rep(jj).stat(1);
        allPostMean(count) = alignStruct(ii).align.rep(jj).stat(3);
    end
end

% Changed selection criteria for noisy stim (alignProtocolDataByTalbe used
% mean difference from median baseline (either pre or total)) 

histFitData = [allPreMean; allPostMean];
fitNames = {'pre'; 'post'};

% fitting 1 and 2 gaussians since sometimes the noise creates a bimodel
% dist (like SPFR for t4new cell 6 amd 7)

modBaseLine = zeros(1,2); 
% this part of the function was designed for protocols with more
% conditions. it should be rewritten
for ii=1:size(histFitData,1)
    
    [histC, histP] = histcounts(histFitData(ii,:), numEdges);
    midP = (histP(1:end-1) + histP(2:end))/2;
    [histFit1, GoF1] = fit(midP', histC', 'gauss1', 'lower', [0.1, min(histFitData(:)), 0.1]);
    if numEdges > 6 
        [histFit2, GoF2] = fit(midP', histC', 'gauss2');
    else
        GoF2.adjrsquare = 0; % so it will skip line 175
    end
    [~, maxHistCI] = max(histC); 
    modBaseLine(ii) = midP(maxHistCI); 
    
    % if there are 2 gaussian fits - low thresh is calculated by lower dist
    % and high by second dist (average is actually more robust this way)  -
    %
    %
    %           !!!      NEED TO CHECK THIS HOLDS  !!!
    %
    %
    % added numEdges > 10 since if there are not enough stim, 2 gaussian
    % estimate is very unreliable
    if GoF1.adjrsquare < 0.8 && GoF2.adjrsquare - 0.2 > GoF1.adjrsquare && numEdges > 10  % to make sure there is a significant difference
        baseThresh(ii).low = histFit2.b1 - sigFac * (histFit2.c1 / sqrt(2)); 
        baseThresh(ii).high = histFit2.b2 + sigFac * (histFit2.c2 / sqrt(2));
        warning(' !!!  fit 2 gaussians %s baseline !!! \n', fitNames{ii})
        baseThresh(ii).adjRSq = GoF2.adjrsquare; 
        if baseThresh(ii).high - baseThresh(ii).low < 2 * sigRng
            error('didnot perform the same correction for 2 hist as I did below')
        end
        
    else
        if sigFac * (histFit1.c1 / sqrt(2)) > sigRng % will allow for a 9mV range
            baseThresh(ii).low = histFit1.b1 - sigFac * (histFit1.c1 / sqrt(2));
            baseThresh(ii).high = histFit1.b1 + sigFac * (histFit1.c1 / sqrt(2));
            baseThresh(ii).adjRSq = GoF1.adjrsquare;
        else
            baseThresh(ii).low = histFit1.b1 - sigRng;
            baseThresh(ii).high = histFit1.b1 + sigRng;
            baseThresh(ii).adjRSq = GoF1.adjrsquare;
        end
    end
    
end

tDiff = zeros(1,2);
for ii=1:2 
    tDiff(ii) = baseThresh(ii).high - baseThresh(ii).low; 
end

baseLab = {'pre', 'post'};

if nargin == 3 && any(tDiff > 35)
    for ii=1:2 
        % use the baseSDStruct to asses Variance and not mean 
        
        preBT(ii).low = baseThresh(ii).low; 
        preBT(ii).high = baseThresh(ii).high; 
        
        preMean = (preBT(ii).low + preBT(ii).high)/ 2; 
        baseHRange = (baseSDStruct(ii).high - baseSDStruct(ii).low) / 2; 
        
        
        baseThresh(ii).low = (baseThresh(ii).adjRSq * preBT(ii).low + baseSDStruct(ii).adjRSq * (preMean - baseHRange)) ...
                              / (baseThresh(ii).adjRSq + baseSDStruct(ii).adjRSq); 
        fprintf('changed low %s baseline from %4.2f to %4.2f \n', baseLab{ii}, preBT(ii).low, preMean - baseHRange)
        
        baseThresh(ii).high = (baseThresh(ii).adjRSq * preBT(ii).high + baseSDStruct(ii).adjRSq * (preMean + baseHRange)) ...
                              / (baseThresh(ii).adjRSq + baseSDStruct(ii).adjRSq);
        fprintf('changed high %s baseline from %4.2f to %4.2f \n', baseLab{ii}, preBT(ii).high, preMean + baseHRange)
    end
end

%checking for the second time (relevant mainly for grating which have very
%few stim (therefore SD estimate is unreliable) 
for ii=1:2 
    tDiff2(ii) = baseThresh(ii).high - baseThresh(ii).low; 
end

if any(tDiff2 > 35) 
    warning('baseline SD estimate unreliable - changed baseline to baseSDStruct')
    assert(nargin == 3, 'baseline SD estimate unreliable - need baseSDStruct input')
    for ii=1:length(baseSDStruct)
        baseHRange = (baseSDStruct(ii).high - baseSDStruct(ii).low) / 2; 
        baseThresh(ii).low = mean(modBaseLine) - baseHRange;
        baseThresh(ii).high = mean(modBaseLine) + baseHRange;
    end
end

if isfield(pStruct.inputParams, 'numCyc')
    fftThresh = [1.5, 0.3]/2; % for longer stim like gratings and flicker
    postTestFlag = 0; % since there can be an OFF resp
else
    fftThresh = [1.5, 0.3]; % determined empirically
    postTestFlag = 1; 
end

stimToCorrect = [];
%finds repeats to exlucde from mean calculation
for ii=1:length(alignStruct)
    for jj=1:length(alignStruct(ii).align.rep)
        
        if postTestFlag
            if  alignStruct(ii).align.rep(jj).stat(1) > baseThresh(1).high || ...
                alignStruct(ii).align.rep(jj).stat(3) > baseThresh(2).high || ...
                alignStruct(ii).align.rep(jj).stat(1) < baseThresh(1).low % no low for post due to persistant inhibition
                
                alignStruct(ii).exclude = [alignStruct(ii).exclude, jj];
                stimToCorrect = [stimToCorrect, ii];
            end
        else
            if  alignStruct(ii).align.rep(jj).stat(1) > baseThresh(1).high || ...
                alignStruct(ii).align.rep(jj).stat(1) < baseThresh(1).low % no low for post due to persistant inhibition
                
                alignStruct(ii).exclude = [alignStruct(ii).exclude, jj];
                stimToCorrect = [stimToCorrect, ii];
            else
                % fft doesnt deal with it since this noise was very fast
                % and needed to be excluded before fft analysis
                if sum(alignStruct(ii).align.rep(jj).data(:,2) == 0) > 10  % excluding strong noise repeats (happen very rarely in G4)
                    alignStruct(ii).exclude = [alignStruct(ii).exclude, jj];
                    stimToCorrect = [stimToCorrect, ii];
                end
            end    
                
        end
            
        if any(alignStruct(ii).align.rep(jj).statFFT - fftThresh > 0)
            
            alignStruct(ii).exclude = [alignStruct(ii).exclude, jj];
            stimToCorrect = [stimToCorrect, ii];
            warning('Stim %d rep %d has high freq noise', ii,jj)
            
        end
        
    end
end

stimToCorrect = unique(stimToCorrect); %since repeats in the same stim can appear multiple times
stimToCorrInds = [];

for ii=1:length(stimToCorrect)
    
    remRep = alignStruct(stimToCorrect(ii)).exclude;
    relReps = setdiff(1:length(alignStruct(stimToCorrect(ii)).align.rep), remRep);
    
    switch length(relReps)
        
        case 0
            
            warning('!!!        Stim %d was REMOVED due to noise        !!!', stimToCorrect(ii))
            alignStruct(stimToCorrect(ii)).align.mean = [];
            alignStruct(stimToCorrect(ii)).align.medMean = [];
            stimToCorrInds = [stimToCorrInds, stimToCorrect(ii)];
            
        case 1
            
            warning('Stim %d has one valid repeat \n', stimToCorrect(ii))
            
            alignStruct(stimToCorrect(ii)).align.mean = alignStruct(stimToCorrect(ii)).align.rep(relReps).data;
            alignStruct(stimToCorrect(ii)).align.meanPos = alignStruct(stimToCorrect(ii)).align.rep(relReps).pos;
            stimToCorrInds = [stimToCorrInds, stimToCorrect(ii)];
            
        otherwise
            
            if iscell(relVarNames) && length(relVarNames) > 1
                posVal = alignTable{stimToCorrect(ii), relVarNames{1}};
            else
                posVal = alignTable{stimToCorrect(ii), relVarNames};
            end
            
            corrInd = alignTable{stimToCorrect(ii), 'relInd'};
            alignMean = getAlignedStimDataByTableExcludeG4(pStruct, corrInd, posVal, relReps, parSt);
            alignStruct(stimToCorrect(ii)).align.mean = alignMean.mean;
            alignStruct(stimToCorrect(ii)).align.medMean = alignMean.medMean;
            alignStruct(stimToCorrect(ii)).align.spikRas = alignMean.spikRas;
            alignStruct(stimToCorrect(ii)).align.meanPos = alignMean.meanPos;
            alignStruct(stimToCorrect(ii)).numReps = alignStruct(stimToCorrect(ii)).numReps - length(alignStruct(stimToCorrect(ii)).exclude);
            
            fprintf('Stim %d has %d valid repeats \n', stimToCorrect(ii), length(relReps))
            stimToCorrInds = [stimToCorrInds, stimToCorrect(ii)];
    end 

    



end

if nargout == 3
    varargout{1} = stimToCorrInds;
elseif nargout == 4
    varargout{1} = stimToCorrInds;
    varargout{2} = baseThresh;
end



end



