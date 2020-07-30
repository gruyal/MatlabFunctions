function mergedStruct  = mergeProtocolsWOverlapStim(pStruct1, pStruct2, disregardSpanFlag)

% function mergedStruct  = mergeProtocolsWOverlapStim(pStruct1, pStruct2)
%
% This function is designed to merge protocols from same createProtocols
% functions that have some overlapping stimuli. It is a modification of
% mergeProtocolsWUniqueStim.
%
%
% INPUT
% pStruct1/2 -          protocolStruct from the same CreateProtocolFunctions
%                       with some overlaping parameters. Note! stimuli dont
%                       need to be uniquely identified by its gratingTable
%                       properties in the combined structure. 
%                       NOTE! tables should still have the same variable
%                       names (generated from same createProtocols
%                       function)
% disregardSpanFlag-    logical (optional). Whether to disregard span when comparing
%                       tables (since it is sometimes just a reference)
%                       defualt = 1
%
%
% OUTPUT 
% mergedStruct -    merged protocolStruct with a combined gratingTable and
%                   stim whose relInds have been appropriately updated
%
%                   Also adds a field .origProt to each stim so that they could be separated
%                   back to thier original protocols and compared (if there is a major
%                   change)
%
%                   NOTE! if one protocol already has a origProt field (since it was already merged) 
%                   it should be the first input. Otherwise it will report
%                   an error

if nargin < 3
    disregardSpanFlag = 1; 
end

%since in some stim the second position is 1 and in others it is the same
%as the first. This differentiates between the two conditions
allP1Inds = vertcat(pStruct1.stim.relInds);
allP2Inds = vertcat(pStruct2.stim.relInds);

pos2P1Inds = unique(allP1Inds(:, 2));
pos2P2Inds = unique(allP2Inds(:, 2));

pos2P1Flag=0;
pos2P2Flag=0;

if all(pos2P1Inds == 1) %use all in case pos2P1IInds is a vector
    pos2P1Flag = 1; 
end

if all(pos2P2Inds == 1) %use all in case pos2P1IInds is a vector
    pos2P2Flag = 1; 
end

assert(pos2P1Flag == pos2P2Flag, 'protocols have different relInds(:,2) convention - check compatibility')

% with repeats
totStim1 = length(pStruct1.stim);
totStim2 = length(pStruct2.stim);

if ~isfield(pStruct1.stim(1), 'origProt')
    for ii=1:totStim1
        pStruct1.stim(ii).origProt = 1; 
    end
    maxProt = 1; 
else
    maxProt = max([pStruct1.stim.origProt]);
end
    
if isfield(pStruct2.stim(1), 'origProt')
    error('protocol with origProt field should be first input in merge')
else
    for ii=1:totStim2
        pStruct2.stim(ii).origProt = maxProt+1;
    end
end
    
mergedStruct = pStruct1; 


%merging tables
firTab = pStruct1.gratingTable; 
secTab = pStruct2.gratingTable; 

varNames = lower(firTab.Properties.VariableNames); 

if disregardSpanFlag
    
    %used strfind and not strcomp since in singleBar prtoocols span is
    %marked uSpan
    spanCell = cellfun(@(x) strfind(x, 'span'), varNames, 'uniformoutput', 0);
    spanInd = find(cellfun(@(x) ~isempty(x), spanCell));
    assert(length(spanInd) <= 1, 'more than one span containing column found')
    
    tabRelInd = setdiff(1:length(varNames), [1, spanInd]);
else
    tabRelInd = 2:length(varNames); % getting rid of index 
end


compTab1 = firTab(:, tabRelInd);
compTab2 = secTab(:, tabRelInd);

if ismember('posSeq', compTab1.Properties.VariableNames) % since setdiff cant deal with sequences 
    compTab1 = removevars(compTab1, 'posSeq');
    compTab2 = removevars(compTab2, 'posSeq'); % since it is the same protocol (currently only appears in movBarShift
end


% finds the rows in table 2 that do not exist in table 1
[~, diffTab2Ind] = setdiff(compTab2, compTab1, 'rows'); % NOTE!! set diff doesn't work with NaNs in the table
diffTab2Ind = sort(diffTab2Ind); % since setdiff changes the order (based on columns in the table
sameTab2Ind = (setdiff(1:height(secTab), diffTab2Ind))';

fprintf('Found %d same stimuli and %d unique \n', length(sameTab2Ind), length(diffTab2Ind))

% for old and new indices
tab2IndIn1 = zeros(length(sameTab2Ind),2);
relTab2 = compTab2(sameTab2Ind, :); 

for ii=1:length(sameTab2Ind)
    newInd = find(ismember(compTab1, relTab2(ii,:)));
    oldInd = sameTab2Ind(ii);
    tab2IndIn1(ii, :) = [oldInd, newInd];
end
    
uniStimTab = secTab(diffTab2Ind, :);

numStim1 = height(firTab);
numUStim2 = height(uniStimTab);

mergedStruct.gratingTable.oldIndex = firTab.index; 
uniStimTab.oldIndex = uniStimTab.index;
uniStimTab.index = ((1:numUStim2) + numStim1)'; 

mergedStruct.gratingTable = vertcat(mergedStruct.gratingTable, uniStimTab);

for ii=1:totStim2
    
    tempInd = pStruct2.stim(ii).relInds(1);
    
    if ismember(tempInd, uniStimTab.oldIndex)
        
        newTempInd = uniStimTab.index(uniStimTab.oldIndex == tempInd);
        
    elseif ismember(tempInd, tab2IndIn1(:,1))
        
        newTempInd = tab2IndIn1(tab2IndIn1(:,1) == tempInd, 2);
        
    else
        error('cant find index in unique or same stimuli')
    end
    
    % since in some protocols the second position is 1 and in others it is
    % equal to the first
    if pos2P1Flag 
        pStruct2.stim(ii).relInds = [newTempInd, 1, 1, 1];
        
    else
        pStruct2.stim(ii).relInds = [newTempInd, newTempInd, 1, 1];
        
    end
    
end



mergedStruct.stim(totStim1+1:totStim1+totStim2) = pStruct2.stim; 

% remove oldIndex so it will not cause problems when applied several times
mergedStruct.gratingTable.oldIndex = [];

end