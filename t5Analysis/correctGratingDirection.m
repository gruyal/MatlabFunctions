function newCellModRep = correctGratingDirection(cellModReport)

% function newCellModRep = correctGratingDirection(cellModReport)
%
% This function fixes the bug in the code before running model8. for cells
% with a flipped PD the table was wrong, This function loads the
% cellModelReport and fixes when needed

allProtTypes = unique(cellModReport.totTable.protType);
modDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T5recordingSummaryAndAnalysis/modelFiles8/';

newCellModRep = cellModReport;

if ismember(1, allProtTypes)
    
    % loading a iteration file to determine direction 
    tempIter = organizingClusterData(cellModReport.cellNum, 555, modDir); %just random iter
    relDir = tempIter.dir;
    
    if relDir == -1
        tempDir = newCellModRep.totTable.direction(newCellModRep.totTable.protType == 1);
        newCellModRep.totTable.direction(newCellModRep.totTable.protType == 1) = ~tempDir;
    end
    
end



end
    
    