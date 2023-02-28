function align = getAlignedStimDataByTableG4(pStruct, gratingInd, posVals, parStruct)

% function alignStimMat = getAlignedStimDataByTableG4(pStruct, posVals)
%
% This function is a modification of getAlignedStimDataByTable2 which takes
% into account the change in the position fucntion use in G4 code (due to the much
% higher genFreq)
% previous modification added another posVal to mark the end of the stimulus (so that
% calculations of baseline deviations would take the post-stimulus time
% into account also ).
% Also added a blip removal code that replaced electrical blip (changes of
% more than 2.5mV in one sample) with linear interpolation of thier edges
%
% Modifed March 2019 
%
% This function is similar to getAlignedStimData only uses the preexisting gratingTable structure. 
% The original gratingTable needs to be augmented with variables that
% describe the presentation frames of the stim. Data is aligned to given var 
% so that time Zero is when alignVar frame was presented and arranges it
% in a matrix (makes repeats into the same lengths)
%
% INPUT
% pStruct -             protocolStruct w/ stim field and data in it
% gratingInd -          index variable in the gratingTable to be aligned
%                       if a single number function assumes only
%                       gratingTable was used. if 4 number, refers to
%                       stim.relInds directly
% posVals -             Changed in new version to include 3
%                       values alaways. 
%                       1 - posVal to align to 
%                       2 - posVal to calculate preBaseline from 
%                       3 - posVal to calculate postBaseline from (after
%                       timeBuff delay)
%                       1 and second posVals can be the same value
%                       old version description: 
%                       Value for position function that is reported in
%                       data{2} of each stim. Should be the same for all
%                       repeats. If posVal is of length 2 the second
%                       position is used to calculate baseline (baseline
%                       would be defined from first sample to appropriate
%                       sample of second val) - used for protocols like
%                       movingbar where to aligned var in not appear
% ParStruct             (optional) parameter structure with 2 possible
%                       fields:
% .posRefs -            length 3 string (optional). References for which part of consecutive posVal runs to use. 
%                       Since with the new use of position fucntion posVal appears in run.
%                       's' - start of a run (used usually for appearance of stim). 
%                       'm' - middle of a run (used for when the alignment
%                             is to the middle of a stimulus run)
%                       'e' - end of a run (used for the disappearance of a
%                       stimulus). Default is 'sse'. 
%                       Note! do not use "xxx" for specifying string. only 'xxx'
% .medSpkPar -          parameters for the medFiltAndFindSpikes function
%                       (medWin and threshV)
% is one of or both of these fields are not given, defaults will be used
%
%   Note!   alignedVar colum in the table is not generated automaticcaly with the protocol. 
% It should be enetered manually and given in the same units the controller uses (first
% frame is 0)
%
% OUTPUT
% 
% align -               structure with the following fields:
%
% .rep
%   .data -             stim data aligned so posVal timing is zero. 
%                       Also data now is an NX4 matrix with first column
%                       time, second original data, third medfilt data (to
%                       remove spikes) and fourth spike raster 
%   .pos -              posFunc values aligned so posVal timing is zero.
%                       First column is converted from time into index of the data time
%   .stat -             intended for QC. first number is preStimulus mean
%                       and second is entire repeat mean. 
%                       Modified: in this version of the function the third
%                       
%   .statFFT -          results from calcFreqPowerInRelSpectra to eliminate
%                       repeats with high frequency noise. Each number is
%                       the max power within the range of frequencies
%                       given (in fftRanges). 
% .mean -               mean across repeats. To generate a matrix, equal number of samples
%                       are taken before posVal timing, and minimal after. 
% .meanPos -            Position indices after they have been corrected for
%                       the number of smaples dropped from each repeat (to generate the mean)

medWinDef = 351;
threshVDef = 4.5;
smWinDef = 21; 

if nargin < 4
    posRefs = 'sse';
    medSpkPar.medWin = medWinDef;
    medSpkPar.threshV = threshVDef;
else
    if isfield(parStruct, 'posRefs')
        posRefs = parStruct.posRefs;
        assert(length(posRefs) == 3, 'posRefs should have a length of 3')
    else
        posRefs = 'sse';
    end
    
    if isfield(parStruct, 'medSpkPar')
       medSpkPar = parStruct.medSpkPar;
       assert(isfield(medSpkPar, 'medWin'), 'medSpkPar is missing medWin field')
       assert(isfield(medSpkPar, 'threshV'), 'medSpkPar is missing threshV field')
       assert(isfield(medSpkPar, 'smWin'), 'medSpkPar is missing smWin field')
       assert(floor(medSpkPar.medWin) == ceil(medSpkPar.medWin), 'medWin should be a whole number')
       assert(floor(medSpkPar.smWin) == ceil(medSpkPar.smWin), 'medWin should be a whole number')
    else
        medSpkPar.medWin = medWinDef;
        medSpkPar.threshV = threshVDef;
        medSpkPar.smWin = smWinDef;
    end

end




assert(length(posVals) == 3, 'This function requires 3 position values, align, start and end')
assert(posVals(2) <= posVals(3), 'posVal 2 must be smaller than posVal 3')


posVal = posVals(1); 
baseVal = posVals(2); 
endVal = posVals(3); 

for ii=1:3
	assert(ismember(posRefs(ii), {'s', 'm', 'e'}), 'posRefs should be either s, m, or e')
end

posRef = posRefs(1);
baseRef = posRefs(2);
endRef = posRefs(3);

timeBuff = 250; % in ms, added to endTime to calculate baseline after response decayed for timeBuff ms


% relCh = 2; % voltage channel (switched between 2 and 3)
relCh = 3; % voltage channel 

% fprintf('extracting data from channel %d \n', relCh)

datTomV = 10; % factor to multifpy data
timeToms = 10^-3; % converts timing stamps to ms (clock @ 1MHz) 

fftRanges = [58, 62; 3680, 3720]; % empirically found to be 2 ranges for noise (don't know where the second range is coming from


if length(gratingInd) == 1
    indsSt = getStimInds(pStruct, [gratingInd, nan, nan, nan]);
elseif length(gratingInd) == 4
    indsSt = getStimInds(pStruct, gratingInd);
else
    error('gratingInd should be either an index to gratingTable or relInds vector')
end

assert(length(indsSt) == 1, 'values given in getStimIndsInput are not specific to one stimulus')

allPos = pStruct.stim(indsSt(1).inds(1)).data{2}(:,2);


assert(ismember(posVal, allPos), 'posVal is not included in position values for the specified stim')
assert(ismember(baseVal, allPos), 'baseVal is not included in position values for the specified stim')
assert(ismember(endVal, allPos), 'endVal is not included in position values for the specified stim')

tempSt = pStruct.stim(indsSt.inds);
numRepInd = zeros(1,length(tempSt));
for ii=1:length(tempSt)
    numRepInd(ii) = ~isempty(tempSt(ii).data);
end


numReps = find(numRepInd, 1, 'last');

% relPosInds = zeros(1, numReps);
% postPosLen = relPosInds;

count=0;
for ii=1:numReps
   
    tempAll = pStruct.stim(indsSt(1).inds(ii)).data;
    if isempty(tempAll)
        fprintf('getAlignedStimDataByTable: No data for stim index %s repeat %d \n', num2str(gratingInd), ii)
        continue
    end
    
    count=count+1; % to account for empty repeats
    tempPos = double(tempAll{2});
    tempDat = tempAll{1};
    
    if max(tempPos(:,2)) > 30000 % sign that something FUBARed in the controller and therefore cant use the data
        align.rep(count).pos = [];
        align.rep(count).baseTime = 0;
        align.rep(count).endTime = 100;
        align.rep(count).data = [tempDat(:,1), zeros(length(tempDat(:,1)), 3)];
        align.rep(count).stat = zeros(1,3);
        align.rep(count).statFFT = 0;
        
         warning('posVec in stim %s rep %d is corrupt - stimulus rejected', num2str(gratingInd), ii)
        
        if ii==numReps
            tempAll = pStruct.stim(indsSt(1).inds(ii-1)).data;
            tempPos = double(tempAll{2});
        end
        
        continue
    end
    
    
    switch posRef
        case 's'
            posInd = find(tempPos(:,2) == posVal, 1, 'first');
            relTime = tempPos(posInd, 1);
        case 'e'
            % since this is the last frame +1 would move to the beginning of back to background
            posInd = find(tempPos(:,2) == posVal, 1, 'last') + 1; 
            relTime = tempPos(posInd, 1);
        case 'm'
            tempInd1 = find(tempPos(:,2) == posVal, 1, 'first');
            tempInd2 = find(tempPos(:,2) == posVal, 1, 'last');
            posInd = floor((tempInd1+tempInd2)/2);
            relTime = tempPos(posInd, 1);
    end
    
    switch baseRef
        case 's'
            tempInd = find(tempPos(:,2) == baseVal, 1, 'first');
            relBaseTime = tempPos(tempInd, 1);
        case 'e'
            % since this is the last frame +1 would move to the beginning of back to background
            tempInd = find(tempPos(:,2) == baseVal, 1, 'last') + 1; 
            relBaseTime = tempPos(tempInd, 1);
        case 'm'
            tempInd1 = find(tempPos(:,2) == baseVal, 1, 'first');
            tempInd2 = find(tempPos(:,2) == baseVal, 1, 'last');
            tempMidI = floor((tempInd1+tempInd2)/2);
            relBaseTime = tempPos(tempMidI, 1);
    end
            
    switch endRef
        case 's'
            tempInd = find(tempPos(:,2) == endVal, 1, 'first');
            relEndTime = tempPos(tempInd, 1);
        case 'e'
            % since this is the last frame +1 would move to the beginning of back to background
            tempInd = find(tempPos(:,2) == endVal, 1, 'last') + 1; 
            relEndTime = tempPos(tempInd, 1);
        case 'm'
            tempInd1 = find(tempPos(:,2) == endVal, 1, 'first');
            tempInd2 = find(tempPos(:,2) == endVal, 1, 'last');
            tempMidI = floor((tempInd1+tempInd2)/2);
            relEndTime = tempPos(tempMidI, 1);
    end
    

    timeCh = (tempDat(:, 1) - relTime) * timeToms;
    relBTConv = (relBaseTime - relTime) * timeToms;
    relETConv = (relEndTime - relTime) * timeToms;
    dataCh = tempDat(:, relCh) * datTomV;
    posTimeCh = (tempPos(:, 1) - relTime) * timeToms;
    posDat = tempPos(:,2);
    posInDatInd = zeros(length(posTimeCh), 1);
    for jj=1:length(posTimeCh)
        tempInd = find(timeCh - posTimeCh(jj) > 0, 1, 'first');
        if ~isempty(tempInd)
            posInDatInd(jj) = tempInd;
        else
            if jj==length(posTimeCh) % had to error that time is recorded in posDat but missing in data
                posInDatInd = posInDatInd(1:end-1);
                posDat = posDat(1:end-1);
                warning('missing last timestamp in posDat - deleted last position')
            end
        end
    end
    
    relPosInds(count) = posInDatInd(posInd);
    postPosLen(count) = length(dataCh) - relPosInds(count);
    
    align.rep(count).pos = [posInDatInd, posDat];
    align.rep(count).baseTime = relBTConv;
    align.rep(count).endTime = relETConv;
    
    diffDat = diff(dataCh);
    
    % removing electrical blips from the recordings - fast event that are
    % clearly from an outside noise source - replacing it with linear
    % interpolation from the edges of the blip 
    if any(abs(diffDat) > 2.5) % a jump of more than 2.5mV between samples 
        
        sampThresh = 0.5; % in mV - determined empirically
        absDiffVec = smoothdata(abs(diffDat), 'movmean', 11);
        [~, blipPos] = max(absDiffVec);
        preBlipInd = find(absDiffVec(1:blipPos) < sampThresh, 1, 'last');
        postBlipInd = blipPos + find(absDiffVec(blipPos:end) < sampThresh, 1, 'first');
        
        if isempty(preBlipInd)
            preBlipInd = 1;
        end

        if isempty(postBlipInd)
            postBlipInd = length(absDiffVec);
        end
        
        if postBlipInd - preBlipInd > 200
            warning('blip in stim %s rep %d longer than expected - check', num2str(gratingInd), ii)
            dataCh(preBlipInd:postBlipInd) = 0; % to mark it for removal
        else
            corrVal = linspace(dataCh(preBlipInd - 1), dataCh(postBlipInd + 1), postBlipInd - preBlipInd +1);
            dataCh(preBlipInd:postBlipInd) = corrVal;
        end
        
        fprintf('Blip removed from stim %s repeat %d \n', num2str(gratingInd), ii);
        
    end
    
    [medDatCh, spkDatCh] = medFiltAndFindSpikes(dataCh, medSpkPar); 
    align.rep(count).data = [timeCh, dataCh, medDatCh, spkDatCh];
    
    preStimInd = find(timeCh > relBTConv, 1, 'first');
    % to compute baseline after the spikes are removed
    preStimDat = medDatCh(1:preStimInd);
%     gratingInd
    postStimInd = find(timeCh > relETConv + timeBuff, 1, 'first');
    % same as above
    postStimDat = medDatCh(postStimInd:end);
    
    fftRes = calcFreqPowerInRelSpectra(dataCh, fftRanges);
    
    align.rep(count).stat = [mean(preStimDat), mean(medDatCh), mean(postStimDat)];
    align.rep(count).statFFT = fftRes;  
    
end





minPre = min(relPosInds);
minPost = min(postPosLen);

meanData = zeros(minPre+minPost, length(align.rep));
meanTime = meanData;
medMeanData = meanData;
spikeRast = zeros(minPre+minPost, length(align.rep));
meanPosIndx  = nan(size(tempPos, 1), length(align.rep));

for ii=1:length(align.rep) % again, in case there is an empty rep
    
    if isempty(align.rep(ii).pos) % to deal with the new error of posVec with values up to 65K
        continue
    end
    
    posLen(ii) = size(align.rep(ii).pos, 1);
    startIdx = relPosInds(ii)-minPre+1;
    stopIdx = relPosInds(ii) + minPost;
    % changed to 3 so that mean is calculated on the medfilt trace
    meanData(:, ii) = align.rep(ii).data(startIdx:stopIdx, 2);
    medMeanData(:, ii) = align.rep(ii).data(startIdx:stopIdx, 3);
    spikeRast(:, ii) = align.rep(ii).data(startIdx:stopIdx, 4);
    meanTime(:, ii) = align.rep(ii).data(startIdx:stopIdx, 1);
    tempPosAl = align.rep(ii).pos(:,1) - startIdx + 1;
    if length(tempPosAl) ~= size(tempPos,1)
        warning('Not all pos vector have the same length stim %s', num2str(gratingInd))
    end
    meanPosIndx(1:length(tempPosAl),ii) = tempPosAl; % since for mean each repeat is chunked differently
end

posInd = find(size(meanPosIndx,1) == posLen, 1, 'first');
tempPos = align.rep(posInd).pos;

align.mean = [nanmean(meanTime, 2), nanmean(meanData, 2)];
align.medMean = [nanmean(meanTime, 2), nanmean(medMeanData, 2)];
align.spikRas = spikeRast;
align.meanPos = [round(median(meanPosIndx, 2)), tempPos(:,2)]; %since posValues are the same for all reps
align.meanBaseTime = mean([align.rep.baseTime]);

end
