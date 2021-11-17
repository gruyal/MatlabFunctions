function protStruct = loadProtocolFromCellsStructure(cellsStruct, cellNum, protName, protNum)

% function loadProtocolFromCellsStructure(cellNum, protName, protNum)
%
% This function loads the protocolStruct for the specified protocol in the
% given cell. 
%
% INPUT
%
% cellsStruct -         Structure with protocols and file names for each
%                       recorded cell. Generated seperately for T4 and T5
%                       in assigningT5/Cells2Struct scripts
% cellNum -             Integer. Given cell number from cellsStruct
% protName -            (optional) String. name of the protocol from that cell. If not
%                       given function asks for the desired based on
%                       available protocols for that cell. 
% protNum -             (optional) Integer. number of the desired protocol.
%                       in case several are available. If not given the
%                       function suggests available numbers. 
%
% OUTPUT
% 
% protStruct -          loaded protocolStruct


% directory with all the subdirectories for all the cells
baseDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp'; 

numCells = length(cellsStruct);

assert(ismember(cellNum, 1:numCells), 'cell number out of range')

tempCell = cellsStruct(cellNum);
relDir = tempCell.dir; 
fdNames = fieldnames(tempCell);
fdNames = fdNames(2:end); % removes the dir field
relFieldsInds = cellfun(@(x) ~isempty(getfield(tempCell, x)), fdNames);
relNames = fdNames(relFieldsInds);

if nargin < 3
    
    for ii=1:length(relNames)
        fprintf(' %d: %s \n', ii, relNames{ii})
    end
    protInd = input('Please enter the desired protocol by number: \n'); 
    assert(ismember(protInd, 1:length(relNames)), 'protocol number out of range')
    protName = relNames{protInd};
    
else
    
    assert(ismember(protName, relNames), 'protName does not exist for the given cell')
    
end

relProt = getfield(tempCell, protName);
        
if nargin < 4
    
    if length(relProt) == 1
        
        protNum = 1; 
        
    else
      
        for ii=1:length(relProt)
            fprintf('%d: %s \n', ii, relProt(ii).properties)
        end
        protNum = input('Please enter the desired protocol by number: \n'); 
        assert(ismember(protNum, 1:length(relProt)), 'protocol number out of range')  
        
    end
    
else
    
    assert(ismember(protNum, 1:length(relProt)), 'protocol number out of range')  
    
end
            
fileName = fullfile(baseDir, relDir, relProt(protNum).fileName);

fprintf('loading file: %s \n', relProt(protNum).fileName)
    
    
load(fileName)

if contains(relProt(protNum).fileName, 'StructAO')
    protStruct = protocolStructAO;   
else
    protStruct = protocolStruct;   
end
    
    
    
    
end
    
    
