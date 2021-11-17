function [cellsWProtInds varargout] = getCellInd4Protocol(allCellSt, protName)

% function cellsWProtInds = getCellInd4Protocol(allCellSt, protName)
%
% this function get the indicies of all the cells that have data for the
% required protocol
%
%
% INPUT
%
% allCellSt -           structure listing all the protocols for each cell
%                       generated w/ assingningT5Cell2Protocols 
%
% protName -            string. name of the desired protocol


assert(isfield(allCellSt, protName), 'structure does not contain the required prtoocol')


cellsWProtInds = nan(1, length(allCellSt));
count = 0;
outSt = struct;

for ii=1:length(allCellSt)
    
    tempProtInfo = getfield(allCellSt(ii), protName); 

    cellsWProtInds(ii) = ~isempty(tempProtInfo);
    
    if ~isempty(tempProtInfo)
        for jj=1:length(tempProtInfo)
            count = count+1;
            outSt(count).properties = tempProtInfo(jj).properties; 
            outSt(count).cellInd = ii;
        end
    end
    
    
end


if nargout > 1
    varargout{1} = outSt;
end



end
    
