function align = getAlignedStimDataByTable2(pStruct, gratingInd, posVals)

% function alignStimMat = getAlignedStimDataByTable2(pStruct, posVals)
%
% This function is a modification of getAlignedStimDataByTable where I have
% added another posVal to mark the end of the stimulus (so that
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
%
% Note! alignedVar colum in the table is not generated automaticcaly with the protocol. 
% It should be enetered manually and given in the same units the controller uses (first
% frame is 0)
%
% OUTPUT
% 
% align -               structure with the following fields:
%
% .rep
%   .data -             stim data aligned so posVal timing is zero
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


assert(length(posVals) == 3, 'This function requires 3 position values, align, start and end')
assert(posVals(2) < posVals(3), 'posVal 2 must be smaller than posVal 3')


posVal = posVals(1); 
baseVal = posVals(2); 
endVal = posVals(3); 

timeBuff = 500; % in ms, added to endTime to calculate baseline after response decayed for timeBuff ms


relCh = 2; % voltage channel (was 3 in old data)
% relCh = 3; % voltage channel (was 3 in old data)

% fprintf('extracting data from channel %d \n', relCh)

datTomV = 10; % factor to multifpy data
timeToms = 10^-3; % converts timing stamps to ms (clock @ 1MHz) 

fftRanges = [58, 62; 3680, 3720]; % empirically found to be 2 ranges for noise (don't know where the second range is coming from


indsSt = getStimInds(pStruct, [gratingInd, nan, nan, nan]);

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
        fprintf('getAlignedStimDataByTable: No data for stim index %d repeat %d \n', gratingInd, ii)
        continue
    end
    
    count=count+1; % to account for empty repeats
    tempPos = double(tempAll{2});
    tempDat = tempAll{1};
    relTime = tempPos(tempPos(:,2) == posVal, 1);
    relBaseTime = tempPos(tempPos(:,2) == baseVal, 1);
    relEndTime = tempPos(tempPos(:,2) == endVal, 1);
    timeCh = (tempDat(:, 1) - relTime) * timeToms;
    relBTConv = (relBaseTime - relTime) * timeToms;
    relETConv = (relEndTime - relTime) * timeToms;
    dataCh = tempDat(:, relCh) * datTomV;
    posTimeCh = (tempPos(:, 1) - relTime) * timeToms;
    posDat = tempPos(:,2);
    posInDatInd = zeros(length(posTimeCh), 1);
    for jj=1:length(posTimeCh)
        posInDatInd(jj) = find(timeCh - posTimeCh(jj) > 0, 1, 'first');
    end
    
    relPosInds(count) = posInDatInd(posDat == posVal);
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
        
        if postBlipInd - preBlipInd > 200
            warning('blip in stim %d rep %d longer than expected - check', gratingInd, ii)
        end
        
        corrVal = linspace(dataCh(preBlipInd - 1), dataCh(postBlipInd + 1), postBlipInd - preBlipInd +1);
        dataCh(preBlipInd:postBlipInd) = corrVal;
        
        fprintf('Blip removed from stim %d repeat %d \n', gratingInd, ii);
        
    end
    
    align.rep(count).data = [timeCh, dataCh];
    
    preStimInd = find(timeCh > relBTConv, 1, 'first');
    preStimDat = dataCh(1:preStimInd);
    
    postStimInd = find(timeCh > relETConv + timeBuff, 1, 'first');
    postStimDat = dataCh(postStimInd:end);
    
    fftRes = calcFreqPowerInRelSpectra(dataCh, fftRanges);
    
    align.rep(count).stat = [mean(preStimDat), mean(dataCh), mean(postStimDat)];
    align.rep(count).statFFT = fftRes;  
    
end





minPre = min(relPosInds);
minPost = min(postPosLen);

meanData = zeros(minPre+minPost, length(align.rep));
meanTime = meanData;
meanPosIndx  = nan(size(tempPos, 1), length(align.rep));

for ii=1:length(align.rep) % again, in case there is an empty rep
    
    posLen(ii) = size(align.rep(ii).pos, 1);
    startIdx = relPosInds(ii)-minPre+1;
    stopIdx = relPosInds(ii) + minPost;
    meanData(:, ii) = align.rep(ii).data(startIdx:stopIdx, 2);
    meanTime(:, ii) = align.rep(ii).data(startIdx:stopIdx, 1);
    tempPosAl = align.rep(ii).pos(:,1) - startIdx + 1;
    if length(tempPosAl) ~= size(tempPos,1)
        warning('Not all pos vector have the same length stim %d', gratingInd)
    end
    meanPosIndx(1:length(tempPosAl),ii) = tempPosAl; % since for mean each repeat is chunked differently
end

posInd = find(size(meanPosIndx,1) == posLen, 1, 'first');
tempPos = align.rep(posInd).pos;

align.mean = [nanmean(meanTime, 2), nanmean(meanData, 2)];
align.median = [mean(meanTime, 2), median(meanData, 2)];
align.meanPos = [round(median(meanPosIndx, 2)), tempPos(:,2)]; %since posValues are the same for all reps
align.meanBaseTime = mean([align.rep.baseTime]);

end
