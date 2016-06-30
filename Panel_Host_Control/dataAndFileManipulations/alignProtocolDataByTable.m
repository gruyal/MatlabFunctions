function alignStruct = alignProtocolDataByTable(pStruct, relVarName)

% function alignStruct = alignProtocolDataByTable(pStruct, relVarName)
%
% This function uses getAlignedStimDataByTable to align the entire
% protocolStruct and outputs an aligned stucutre. 
%
% INPUT
%
% pStruct -             protocolStruct w/.stim, .data and .gratingTable in it
% relVarName -          variable name in according to which data will be
%                       aligned. 
%       Note!!      Name should come from gratingTable variable and that variable should have 
%                   a relevant controller frame number (starting from 0) accrding to which data would be aligned
%
% OUTPUT
%
% alignStruct -         structure with the same fields as output of
%                       getAlignedDataByTable plus the relevant table row

assert(isfield(pStruct, 'gratingTable'), 'pStruct is missing gratingTable field')

stimTable = pStruct.gratingTable;

assert(ismember(relVarName, stimTable.Properties.VariableNames), 'gratingTable does not contain relVarName')
assert(ismember('index', stimTable.Properties.VariableNames), 'gratingTable does not contain index variable')

alignStruct = struct;

for ii=1:height(stimTable)
    
    relTable = stimTable(ii, :);
    relVal = relTable{:, relVarName};
    
    alignStruct(ii).align = getAlignedStimDataByTable(pStruct, relTable.index, relVal);
    alignStruct(ii).table = relTable;
    
end




end

