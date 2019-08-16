function inputDat = extractADCDataFromTDMSStruct(tdmsStruct)

% function xPos = extractPositionFromTDMSStruct(tdmsStruct)
%
% This function exctracts the x position from the TDMS command data. 
% Strongly depends on the TDMS structure and focuses only on dump frame
% commands. 
% Function can be modified to take time stams from all channels (if
% smapling rate is different)
%
% INPUT
% tdmsStruct -  generated by reading log file with TDMS_readTDMSFile
%
% OUTPUT
% inputDat -    recorded voltage and time stamp (only from one channel) from all the ADC# channels


relGroups  = cellfun(@(x) strfind(x, 'ADC'), tdmsStruct.groupNames, 'uniformoutput', 0);
groupInd = find(cellfun(@(x) ~isempty(x), relGroups));

% Making sure ADC channels are read in in the correct order (by name)
adcGroupNames = vertcat(tdmsStruct.groupNames{groupInd});
[~, ord] = sort(adcGroupNames(:,end));
groupInd = groupInd(ord);

% assumes equal sampling rate on all ADC channels 
numSamp = max(tdmsStruct.numberDataPointsRaw);
inputDat = zeros(length(groupInd)+1, numSamp);

for ii=1:length(groupInd)
    
    if ii==1
        timeInd = find(cellfun(@(x) strcmp(x, 'Time'), tdmsStruct.chanNames{groupInd(ii)}));
        timeStInd = tdmsStruct.chanIndices{groupInd(ii)}(timeInd);
        timeStamp = tdmsStruct.data{timeStInd};
        inputDat(ii, :) = timeStamp;
    end
    
    chanInd = find(cellfun(@(x) strcmp(x, 'Voltage'), tdmsStruct.chanNames{groupInd(ii)}));
    dataInd = tdmsStruct.chanIndices{groupInd(ii)}(chanInd);
    vData = tdmsStruct.data{dataInd};
    inputDat(ii+1, :) = vData;
    

end
    


end