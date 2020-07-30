function varargout = plotAlignFlightCLNew(expStruct)

% function varargout = plotAlignFlightCL(expStruct)
% 
%  this funciton plots close loop flight experiments from the same fly
%
% currently alignment marks CL experiments with continuous flight with
% usefulTag
%
% New version of this function has several fixation patterns and uses
% expTable to read pattern properties 

clf

numPos = expStruct.dataCL(1,1).table.numFrames; % in fixation pattern (assmues all fixation patterns have the same number of frames



clData = expStruct.dataCL; 
datSiz = size(clData);
cutOffVal = floor(datSiz(2)/2);

appGroups = cell(datSiz(1), 2); % fixation patterns and begining and end of exp


% useTag = [clData.usefulTag];


for ii=1:datSiz(1)
    
    for jj=1:datSiz(2)
        
        
        
        if ~clData(ii,jj).usefulTag || ~clData(ii,jj).table.presented
            continue 
        end
    
    
        relPos = clData(ii,jj).position(:,2);
        
        % to fold histogram over (since pattern flips in front of fly)
        rotInds = relPos <= floor(numPos/2);
        relPos(rotInds) = relPos(rotInds) + numPos;
        if jj > cutOffVal 
            grpInd = 2; 
        else
            grpInd = 1;
        end

        appGroups{ii, grpInd} = vertcat(appGroups{ii,grpInd}, relPos);
        
    end
    
end

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -0.01, 0.02, datSiz(1));

axh = gobjects(datSiz(1), 1);

plotBins = (1:2:numPos) - 1.5 + floor(numPos/2); 

for ii=1:datSiz(1)
        
    axh(ii) = axes('position', posCell{ii});
    
    histogram(appGroups{ii, 1}, plotBins)
    hold on 
    histogram(appGroups{ii, 2}, plotBins)
    yyLim = axh(ii).YLim; 
    
    line([numPos, numPos], yyLim, 'color', 'r', 'linewidth', 2)
    axh(ii).XLim = floor(numPos/2) + [1, numPos];
    hold off
        
end


if nargout > 0
    varargout{1} = axh;
end
    

end



    
    
    
    
