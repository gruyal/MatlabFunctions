function newExpTable = addRepNumerToExpTable(expTable)

% function newExpTable = addRepNumerToExpTable(expTable)
%
% This function takes the expTable and adds a column for the repeat number
% for each stimulus (defined by pattern ID) 

numPats = unique(expTable.patNum);

newExpTable = expTable; 

newExpTable.repNum = zeros(height(expTable), 1); 

for ii=1:length(numPats)
    
    tempInds = find(expTable.patNum == numPats(ii)); 
    
    newExpTable.repNum(tempInds) = 1:length(tempInds);
    
end
    
    





end

