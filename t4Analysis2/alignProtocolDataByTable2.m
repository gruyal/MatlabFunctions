function [alignStruct, varargout] = alignProtocolDataByTable2(pStruct, relVarNames, baseSDStruct)

% function alignStruct = alignProtocolDataByTable2(pStruct, relVarName)
%
% This function uses getAlignedStimDataByTable to align the entire
% protocolStruct and outputs an aligned stucutre. 
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
% stimToCorrInds -      indices for stim that were deemed inadmissable


if nargin < 3
    for ii=1:2
        baseSDStruct(ii).low = -65;
        baseSDStruct(ii).high = -65;
        baseSDStruct(ii).adjRSq = 0; % weight for this would be zero
    end
end


sigFac = 3.5; % factor to multifly sigma and exclude everything above mu + sigFac * sigma (from pre and post baseline calculation)

assert(isfield(pStruct, 'gratingTable'), 'pStruct is missing gratingTable field')

stimTable = pStruct.gratingTable;

assert(all(ismember(relVarNames, stimTable.Properties.VariableNames)), 'gratingTable does not contain one or all of relVarNames')
assert(length(relVarNames) == 3, 'relVarNames must be of length 3 - align, start and end names (align and start can be the same')
assert(ismember('index', stimTable.Properties.VariableNames), 'gratingTable does not contain index variable')

alignStruct = struct;

for ii=1:height(stimTable)
    
    relTable = stimTable(ii, :);
    relVal = relTable{:, relVarNames};
    
    % need to check for minMot protocols !!
    if any(relVal == -1) && strcmp(relVarNames, 'fAppear') % distinction for minMot protocols with First Bar being different color from second when they are presented in same position
        relVal = relTable{:, 'sAppear'};
        fprintf('Changed relVarName to sAppear for stim %d \n', relTable.index)
    end

    alignStruct(ii).align = getAlignedStimDataByTable2(pStruct, relTable.index, relVal);
    alignStruct(ii).table = relTable;
    alignStruct(ii).exclude = [];
    alignStruct(ii).numReps = length(alignStruct(ii).align.rep);
    
end

totReps = sum([alignStruct.numReps]); 

numEdges = floor(totReps/20 / 10) * 10; % worked for SB that has a lot of repeats 

if numEdges < 12 % since above calculation is for protocols with totReps > 500
    numEdges = floor(totReps/10); 
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

for ii=1:size(histFitData,1)
    
    [histC, histP] = histcounts(histFitData(ii,:), numEdges);
    midP = (histP(1:end-1) + histP(2:end))/2;
    [histFit1, GoF1] = fit(midP', histC', 'gauss1');
    [histFit2, GoF2] = fit(midP', histC', 'gauss2');
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
    if GoF2.adjrsquare - 0.15 > GoF1.adjrsquare && numEdges > 10  % to make sure there is a significant difference
        baseThresh(ii).low = histFit2.b1 - sigFac * (histFit2.c1 / sqrt(2)); 
        baseThresh(ii).high = histFit2.b2 + sigFac * (histFit2.c2 / sqrt(2));
        warning(' !!!  fit 2 gaussians %s baseline !!! \n', fitNames{ii})
        baseThresh(ii).adjRSq = GoF2.adjrsquare; 
        
    else
        baseThresh(ii).low = histFit1.b1 - sigFac * (histFit1.c1 / sqrt(2));
        baseThresh(ii).high = histFit1.b1 + sigFac * (histFit1.c1 / sqrt(2));
        baseThresh(ii).adjRSq = GoF1.adjrsquare; 
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
            stimToCorrInds = [stimToCorrInds, stimToCorrect(ii)];
            
        case 1
            
            warning('Stim %d has one valid repeat \n', stimToCorrect(ii))
            
            alignStruct(stimToCorrect(ii)).align.mean = alignStruct(stimToCorrect(ii)).align.rep(relReps).data;
            alignStruct(stimToCorrect(ii)).align.meanPos = alignStruct(stimToCorrect(ii)).align.rep(relReps).pos;
            stimToCorrInds = [stimToCorrInds, stimToCorrect(ii)];
            
        otherwise
            
            if iscell(relVarNames) && length(relVarNames) > 1
                posVal = stimTable{stimToCorrect(ii), relVarNames{1}};
            else
                posVal = stimTable{stimToCorrect(ii), relVarNames};
            end
                
            alignMean = getAlignedStimDataByTableExclude(pStruct, stimTable{stimToCorrect(ii), 'index'}, posVal, relReps);
            alignStruct(stimToCorrect(ii)).align.mean = alignMean.mean;
            alignStruct(stimToCorrect(ii)).align.meanPos = alignMean.meanPos;
            
            fprintf('Stim %d has %d valid repeats \n', stimToCorrect(ii), length(relReps))
            stimToCorrInds = [stimToCorrInds, stimToCorrect(ii)];
    end 

    



end

if nargout == 2
    varargout{1} = stimToCorrInds;
elseif nargout == 3
    varargout{1} = stimToCorrInds;
    varargout{2} = baseThresh;
end



end


