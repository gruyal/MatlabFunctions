function varargout = plotAlignedResponsesWRepeats(alignStruct, xxRange)

% function varargout = plotAlignedResponsesWRepeats(mmResultStruct)
%
% This function simply plots all the different resposnes with all their
% repeats after they have been aligned. Repeats that have been excluded are
% labelled with a different color
%
% input should be output from alignProtocolDataByTable2
% xxRange is optional 

if nargin < 2
    xxRange = 0;
end


close all

allStim = length(alignStruct); 

numStimPerPlot = ceil(sqrt(allStim)); 
numHigh = 8; 

if numStimPerPlot > numHigh
    numStimPerPlot = numHigh;
end

posCell = generatePositionCell(0.05, 0.975, 0.025, 0.975, 0.01, 0.01, [numStimPerPlot, numStimPerPlot]);
count = 0;

pCol = cbrewer('qual', 'Set1', 4);
repCol = [1,1,1]*0.8;
axh = gobjects(1, allStim);
yyMin = 0;
yyMax = -100;


figure('position', [605, 150, 1400, 1000])

newFigTresh = numStimPerPlot^2;

for ii=1:allStim
    
    if ii > newFigTresh 
        figure('position', [605, 150, 1400, 1000])
        count = 0;
        newFigTresh = newFigTresh + numStimPerPlot^2; 
    end
    
    
    relData = alignStruct(ii).align; 
    
    meanD = relData.mean; 
    count=count+1;
    axh(ii) = axes('position', posCell{count});
    hold on 
    
    for jj=1:length(relData.rep)
        repD = relData.rep(jj).data; 
        
        if ismember(jj, alignStruct(ii).exclude)
            relRepCol = pCol(2,:);
        else
            relRepCol = repCol; 
        end
        
        plot(repD(:,1), repD(:,2), 'linewidth', 2, 'color', relRepCol)
        
        if min(repD(:,2)) < yyMin
            yyMin = min(repD(:,2)); 
        end
        
        if max(repD(:,2)) > yyMax
            yyMax = max(repD(:,2)); 
        end
        
    end

    
    plot(meanD(:,1), meanD(:,2), 'linewidth', 3, 'color', pCol(1,:))
    title(num2str(ii))
    if xxRange == 0
        axh(ii).XLim = [meanD(1,1), meanD(end,1)]; 
    else
        axh(ii).XLim = xxRange;
    end
    axh(ii).XColor = 'none';
    axh(ii).YColor = 'none';
    
    if rem(count, numStimPerPlot) ==1
        axh(ii).YColor = 'k';
    end

    
end

cyMax = ceil(yyMax/10)*10;
cyMin = round(yyMin/10)*10;
yyRange = cyMax - cyMin; 
yyLim = [cyMin - yyRange/10, cyMax + yyRange/10];

for ii=1:allStim
    axh(ii).YLim = yyLim; 
end



if nargout == 1
    varargout{1} = axh;
end






end