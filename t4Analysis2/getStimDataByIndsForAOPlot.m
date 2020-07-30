function [datMat, timMat, aodMat, aotMat] = getStimDataByIndsForAOPlot(pStruct, inds, relCh)

% function [datMat, timMat] = getStimDataByIndsForAOPlot(pStruct, inds)
% 
% This is a modification of getStimDataByInds to include the AO function as
% output (using the same time vector)
%
% This function extracts the relevant data from pStruct (accroding to the
% inds) and generates data and timestamp matrices. 
% Function does not check inputs since it is used internally after checks
% have been made
%
% INPUT
% pStruct -         regular protocol structure with .stim.data fields
% inds -            1XN vector of requested indices
% relCh -           optional. channel to get data from (2 for T5 3 for T4
%                   <older data>. if not given default is 2
%
% OUTPUT
% datMat/timMat -   data and timestamp matrices padded with nans to allow
%                   for mean calculations


if nargin < 3
    relCh = 2; %current channel
end

datToMvConv = 10;
datToMsConv = 10^-3;

tempData = cell(1, length(inds));
tempTimS = cell(1, length(inds));
tempAOD = cell(1, length(inds));
tempAOT = cell(1, length(inds));

for ii=1:length(inds)
    if isempty(pStruct.stim(inds(ii)).data)
        warning('Data field for stim Index %d is empty', inds(ii))
    else
        currD = pStruct.stim(inds(ii)).data{1};
        AOD = pStruct.stim(inds(ii)).data{3};
        tempData{ii} = currD(:, relCh) * datToMvConv;
        tempAOD{ii} = AOD(:,2); 
        tempTimS{ii} = (currD(:, 1) - currD(1,1)) * datToMsConv; % zeros the first timestamp for each trial
        tempAOT{ii} = (AOD(:,1) - currD(1,1)) * datToMsConv;
    end
end

allLen = cellfun(@length, tempData);
allLenAO = cellfun(@length, tempAOD);

datMat = nan(max(allLen), length(inds));
timMat = datMat;

aodMat = nan(max(allLenAO), length(inds));
aotMat = aodMat;

for ii=1:length(allLen)
    datMat(1:allLen(ii), ii) = tempData{ii};
    timMat(1:allLen(ii), ii) = tempTimS{ii};
end

for ii=1:length(allLenAO)
    aodMat(1:allLenAO(ii), ii) = tempAOD{ii};
    aotMat(1:allLenAO(ii), ii) = tempAOT{ii};
end


% if an index refers to empty data (since experiment was aborted); that
% index is excluded
notEmptyReps = ~isnan(nanmean(datMat));

datMat = datMat(:,notEmptyReps);
timMat = timMat(:,notEmptyReps);
aodMat = aodMat(:,notEmptyReps);
aotMat = aotMat(:,notEmptyReps);


end