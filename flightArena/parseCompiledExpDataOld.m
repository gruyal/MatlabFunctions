function parseCompDatSt = parseCompiledExpDataOld(compDat, expOrder)

% function parseCompiledExpData(comDat, expOrder)
%
% This function uses the output of compileExpData together with the
% exp_order file to create a parsed structure of responses by stimulus. 


freqCh = 5; 
expOrder = expOrder'; % should be Stim X Reps (otherwise order is fucked)

relTab = compDat.commTable;

% binning stimuli
stimStInd = find(strcmpi(relTab.commType, 'start-display'));
% stimEndInd = find(strcmpi(relTab.commType, 'stop-display')); % unnecessary since always 1 bigger than start and has one extra
preStimTypInd = find(strcmpi(relTab.commType, 'Set Control Mode'));

stimTypeCell = relTab.commVal(preStimTypInd);
stimTypeVal = cellfun(@str2num, stimTypeCell);

assert(length(stimStInd) == length(stimTypeVal), ...
       'control mode and start display do not have same number of commands')

assert(2*numel(expOrder) == length(stimStInd), ...
       'stimuli number from expOrder does not match start_dispay inds') % times 2 since fixation is interleaved in each trial

timeBuff = 0; % (in secs) in case want to add time before and/or after recorded timestamp    

for ii=1:length(stimStInd)
    
    tempStartTime = double(relTab.commTime(stimStInd(ii)));
    tempEndTime = double(relTab.commTime(stimStInd(ii) + 1));
    
    datStInd = find(compDat.chData(:,1) - tempStartTime - timeBuff > 0, 1, 'first');
    datEndInd = find(compDat.chData(:,1) - tempEndTime + timeBuff > 0, 1, 'first');

    relDat = compDat.chData(datStInd:datEndInd, :);
    
    posStInd = find(compDat.frameData(:,1) - tempStartTime > 0, 1, 'first');
    posEndInd = find(compDat.frameData(:,1) - tempEndTime > 0, 1, 'first');

    posDat = compDat.frameData(posStInd:posEndInd, :);
    
    
    if stimTypeVal(ii) == 1 % open loop exp
        
        tempStimInd = ceil(ii/2); % to avoid counting fixation
        stimID = expOrder(tempStimInd);
        stimOrd = find(expOrder == stimID);
        repID = find(stimOrd == tempStimInd);
        
        parseCompDatSt.dataOL(stimID - 1, repID).data = relDat; % -1 since 1 is fixation 
        parseCompDatSt.dataOL(stimID - 1, repID).position = posDat;
        parseCompDatSt.dataOL(stimID - 1, repID).time = [tempStartTime, tempEndTime]; % based on start and stop display 
        
    else
        
        tempStimInd = ii/2; % should always be round number
        parseCompDatSt.dataCL(tempStimInd).data = relDat; 
        
        parseCompDatSt.dataCL(tempStimInd).position = posDat;
        
    end


end


end


