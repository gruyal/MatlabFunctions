function parseCompDatSt = parseCompiledExpData(compDat, expTable)

% function parseCompiledExpData(comDat, expTable)
%
% This function uses the output of compileExpData together with the
% expTable file to create a parsed structure of responses by stimulus. 
%
% Sameas old version only uses expTable (created when running the
% experiment) instead of exp_order


relTab = compDat.commTable;

modCell = expTable.pType;
modCell = modCell(logical(expTable.presented));
modVec = zeros(length(modCell), 1);

for ii=1:length(modCell)
    if strcmp(modCell{ii}, 'fix')
        modVec(ii) = 4; 
    elseif strcmp(modCell{ii}, 'opto')
        modVec(ii) = 1; 
    end
end

% binning stimuli
stimStInd = find(strcmpi(relTab.commType, 'start-display'));
% stimEndInd = find(strcmpi(relTab.commType, 'stop-display')); % unnecessary since always 1 bigger than start and has one extra
preStimTypInd = find(strcmpi(relTab.commType, 'Set Control Mode'));

% since experiment was marked as started but not as presented
if sum(expTable.presented) < height(expTable) % if experiment was aborted
    stimStInd = stimStInd(1:end-1);
    preStimTypInd = preStimTypInd(1:end-1);
end

stimTypeCell = relTab.commVal(preStimTypInd);
stimTypeVal = cellfun(@str2num, stimTypeCell);

assert(length(stimStInd) == length(stimTypeVal), ...
       'control mode and start display do not have same number of commands')

% assert(height(expTable) == length(stimStInd), ...
%        'stimuli number from expOrder does not match start_dispay inds') % times 2 since fixation is interleaved in each trial
   
assert(all(eq(stimTypeVal, modVec)), 'experiment mode from table and does not match recorded')

timeBuff = 0; % (in secs) in case want to add time before and/or after recorded timestamp    

firstOptoNum = min(expTable.patNum(strcmp(expTable.pType, 'opto')));
firstFixNum = min(expTable.patNum(strcmp(expTable.pType, 'fix')));

for ii=1:height(expTable)
    
    stimTable = expTable(ii,:);
    
    % in case stim wasn't even presented (exp aborted)
    if ~stimTable.presented
        continue
    end
    
    tempStartTime = double(relTab.commTime(stimStInd(ii)));
    tempEndTime = double(relTab.commTime(stimStInd(ii) + 1));
    
    datStInd = find(compDat.chData(:,1) - tempStartTime - timeBuff > 0, 1, 'first');
    datEndInd = find(compDat.chData(:,1) - tempEndTime + timeBuff > 0, 1, 'first');

    relDat = compDat.chData(datStInd:datEndInd, :);
    
    posStInd = find(compDat.frameData(:,1) - tempStartTime > 0, 1, 'first');
    posEndInd = find(compDat.frameData(:,1) - tempEndTime > 0, 1, 'first');
    
    % needed to add this due to lack of delay at the last stim (final stop
    % command is recorded after the last frame is recorded
    if isempty(posEndInd)
        posEndInd = length(compDat.frameData(:,1));
        assert(compDat.frameData(posEndInd, 1) - tempEndTime < 0.1, 'difference too big')
    end

    posDat = compDat.frameData(posStInd:posEndInd, :);
    
    stimID = stimTable.patNum; 
    repID = stimTable.repNum; 
    
    if strcmp(stimTable.pType, 'opto') % open loop exp
        
        structInd = stimID - firstOptoNum + 1; % since fixation patterns are first
        
        parseCompDatSt.dataOL(structInd , repID).table = stimTable; 
        parseCompDatSt.dataOL(structInd , repID).data = relDat; 
        parseCompDatSt.dataOL(structInd , repID).position = posDat;
        parseCompDatSt.dataOL(structInd , repID).time = [tempStartTime, tempEndTime]; % based on start and stop display 
        
    elseif strcmp(stimTable.pType, 'fix') 
        
        structInd = stimID - firstFixNum + 1; % since fixation patterns are first
        
        parseCompDatSt.dataCL(structInd , repID).table = stimTable; 
        parseCompDatSt.dataCL(structInd , repID).data = relDat; 
        parseCompDatSt.dataCL(structInd , repID).position = posDat;
        parseCompDatSt.dataCL(structInd , repID).time = [tempStartTime, tempEndTime]; % based on start and stop display 
        
    end


end


end


