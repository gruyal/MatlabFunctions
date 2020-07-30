function varargout = plotAlignFlightCL(expStruct)

% function varargout = plotAlignFlightCL(expStruct)
% 
%  this funciton plots close loop flight experiments from the same fly
%
% currently alignment marks CL experiments with continuous flight with
% usefulTag

clf

numPos = 192; % in fixation pattern
numGroup = 4; 

clData = expStruct.dataCL; 
useTag = [clData.usefulTag];

groupSiz = ceil(length(clData)/numGroup); % divide into 4 groups to plot throughout experiment
appGroups = cell(1, numGroup);

for ii=1:length(clData)
    
    if ~useTag(ii)
        continue 
    end
    
    grpInd = floor(ii/(groupSiz+0.1)) + 1; 
    relPos = clData(ii).position(:,2);
    
    % to fold histogram over (since pattern flips in front of fly)
    rotInds = relPos <= floor(numPos/2);
    relPos(rotInds) = relPos(rotInds) + numPos;
    appGroups{grpInd} = vertcat(appGroups{grpInd}, relPos);
    
end

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -1, 0.01, numGroup);

axh = gobjects(1,numGroup+1);

plotBins = (1:2:numPos) - 1.5 + floor(numPos/2); 

for ii=1:numGroup
    
    axh(ii) = axes('position', posCell{ii});
        
    histogram(appGroups{ii}, plotBins)
    yyLim = axh(ii).YLim; 
    
    line([numPos, numPos], yyLim, 'color', 'r', 'linewidth', 2)
    axh(ii).XLim = floor(numPos/2) + [1, numPos];

end


if nargout > 0
    varargout{1} = axh;
end
    

end



    
    
    
    
