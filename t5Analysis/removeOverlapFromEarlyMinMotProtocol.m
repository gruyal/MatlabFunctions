function newMMProt = removeOverlapFromEarlyMinMotProtocol(mmProt)

% function newMMProt = removeOverlapFromEarlyMinMotProtocol(mmProt)
%
% this function is designed to remove minMot stimuli in which the 2 bars
% are overlapping. There was a bug in the code before March 2018 that
% casued the first bar to be presented partially occluded. 
% Since March 2018 the code has been fixed and first bar is presented
% fully. 

newMMProt = mmProt; 

ov = mmProt.gratingTable.overlap;

grtInds = mmProt.gratingTable.index(ov ~=1);

allInds = vertcat(mmProt.stim.relInds);

selStimInds = ismember(allInds(:,1), grtInds);

newMMProt.gratingTable = mmProt.gratingTable(grtInds, :); 
newMMProt.stim = mmProt.stim(selStimInds); 

stimRemoved = length(ov) - length(grtInds); 

fprintf('Removed %d stimuli from a total of %d \n', stimRemoved, length(ov))

% reindexing the stim so that they wont miss numbers
if stimRemoved
    
    prevIndex = newMMProt.gratingTable.index; 
    newIndex = (1:length(prevIndex))';
    newMMProt.gratingTable.index = newIndex; 
    
    for ii=1:length(newMMProt.stim)
        
        tempPrevInd = newMMProt.stim(ii).relInds(1); 
        tempNewInd = find(prevIndex == tempPrevInd);
        newMMProt.stim(ii).relInds(1) = tempNewInd; 
        
    end
    
end
    




end
    
