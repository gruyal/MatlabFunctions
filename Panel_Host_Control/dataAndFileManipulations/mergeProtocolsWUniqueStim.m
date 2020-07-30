function mergedStruct  = mergeProtocolsWUniqueStim(pStruct1, pStruct2)

% function mergedStruct  = mergeProtocolsWUniqueStim(pStruct1, pStruct2)
%
% This function is designed to merge protocols from same createXXXProtocols
% functions that have different parameters, such that each stimulus in both
% protocols can still be uniquely identified in the combined gratingTable
%
%
% INPUT
% pStruct1/2 -      protocolStruct from the same CreateProtocolFunctions
%                   with different parameters. Note! each stimulus still
%                   needs to be uniquely identified by its gratingTable
%                   properties in the combined structure
%
%
% OUTPUT 
% mergedStruct -    merged protocolStruct with a combined gratingTable and
%                   stim whose relInds have been appropriately updated
%
%       Note! function does not check if there is overlap. It will simply
%       generate non-unique table enteries
%
mergedStruct = pStruct1; 

% with repeats
totStim1 = length(pStruct1.stim);
totStim2 = length(pStruct2.stim);

%merging tables
numStim1 = height(pStruct1.gratingTable);
numStim2 = height(pStruct2.gratingTable);
secTab = pStruct2.gratingTable; 
secTab.index = secTab.index + numStim1; 

mergedStruct.gratingTable = vertcat(mergedStruct.gratingTable, secTab);

mergedStruct.stim(totStim1+1:totStim1+totStim2) = pStruct2.stim; 

for ii=1:totStim2
    
    mergedStruct.stim(totStim1+ii).relInds = mergedStruct.stim(totStim1+ii).relInds + [numStim1, numStim1, 0, 0];
    
end
    







end