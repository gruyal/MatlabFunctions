function varargout = plotMinMotResultsByTable2Axes(mmResStruct, parNameAndVal, plotParSt)

% function varargout = plotStimStructByTable(mmResStruct, figDim, posDim, colDim, plotSt)
%
% This function is a modification of plotStimStructByTable2Axes, designed specifically
% for output from calcMinMotExtLinCompDiffWandVwTable. 
%
%
% This functions plots the data in mmResStruct based on the separation in the
% different dimensions. mmResStruct should have a .stim field, and should work
% with getStimInds (have either a relInds or a combInds field). 
% 
% INPUTS
% mmResStruct -     protocol structure with .stim field and .data in it
% parNameAndVal -   structure. Should contain an include fields with 4 fields each with 'name' and
%                   'values' subfields.
% .include
%   1 -             (figure) variable name in the mmTable that will be plotted as
%                   spearate figure, and the values to use when plotting
%   2/3 -           (axes 1 & 2) same for axes in the horizontal/vertical direction
%   4 -             (color) same for color.
% .exclude
%                   Also can contain an exclude field (with N fields each
%                   with name and values subfields) to further limit the
%                   data by exclluding the given combinations from the
%                   table
%
% if values is empty it takes all values
%
% plotParSt -       (optional) structure with different parameters that
%                   will affect data plotting 
% .figPar -         fields within .figPar will be applied on each figure.
%                   if only one is given same paraemters will be applied on
%                   all figures generated. 
%                   figPar can contain the following fields:
%   .axesNum -      2 element vector to describe number of axes to plot on
%                   the X and Y dimensitons of each figure. default is
%                   [floor(sqrt(numAxesDim)), ceil(sqrt(numAxesDim))]
%   .axesOrd -      1XnumAxes vector describing the order in which axes are
%                   to be populated. Default is 1:numAxes
%   .cols -         2*numColDimX3 colors to be used in plotting the data
%                   for the figure. default pallete is cbrewer(Paired), to
%                   match mean and individual repeats. If input is just
%                   numColDimX3 mean and indiv are at the same color
%
% if figPar is left empty for a particular figure it will be generated by
% default parameters.
%
% fig/axe/colName can also be none if it not needed


%% verifying input correctness

close all

assert(isfield(mmResStruct, 'mmTable'), 'mmResStruct is missing mmTable field')

relTable = mmResStruct.mmTable;
relTable.none = ones(height(relTable),1);
colTableNames = relTable.Properties.VariableNames;

assert(length(parNameAndVal.include) == 4, 'parNameAndVal structure should have 4 fields (fig, axe1,axe2, and color)')

figName = parNameAndVal.include(1).name;
figVals = parNameAndVal.include(1).values;

axeNames{1} = parNameAndVal.include(2).name;
axeNames{2} = parNameAndVal.include(3).name;

axeVals1 = parNameAndVal.include(2).values;
axeVals2 = parNameAndVal.include(3).values;

colName = parNameAndVal.include(4).name;
colVals = parNameAndVal.include(4).values;

assert(all(ismember([{figName}, axeNames, {colName}], colTableNames)), 'Given name does not match variable names in table')

if ~isempty(figVals)
    relTable = relTable(ismember(relTable{:, figName}, figVals), :);
end

if ~isempty(axeVals1)
    relTable = relTable(ismember(relTable{:, axeNames{1}}, axeVals1), :);
end

if ~isempty(axeVals2)
    relTable = relTable(ismember(relTable{:, axeNames{2}}, axeVals2), :);
end

if ~isempty(colVals)
    relTable = relTable(ismember(relTable{:, colName}, colVals), :);
end

if isfield(parNameAndVal, 'exclude')
    
    for ii=1:length(parNameAndVal.exclude)
        exName = parNameAndVal.exclude(ii).name; 
        exVals = parNameAndVal.exclude(ii).values;
        
        relTable = relTable(~ismember(relTable{:, exName}, exVals), :);
        
    end
    
end

allAx1Vals = relTable{:,axeNames{1}}; 
allAx2Vals = relTable{:,axeNames{2}}; 


figIndVals = unique(relTable{:, figName});
numFigs = length(figIndVals);


for ff=1:numFigs
    tempTab = relTable(relTable{:, figName} == figIndVals(ff), :);
    axeRowIndVals = unique(tempTab{:, axeNames{1}});
    axeColIndVals = unique(tempTab{:, axeNames{2}});
    colIndVals = unique(tempTab{:, colName});
    plotSt.figPar(ff).relAxeRowVal = axeRowIndVals;
    plotSt.figPar(ff).relAxeColVal = axeColIndVals;
    plotSt.figPar(ff).relColVal = colIndVals;
    plotSt.figPar(ff).relNumAxeRow = length(axeRowIndVals);
    plotSt.figPar(ff).relNumAxeCol = length(axeColIndVals);
    plotSt.figPar(ff).relNumCol = length(colIndVals);
    plotSt.figPar(ff).axesNum = [length(axeColIndVals), length(axeRowIndVals)];
    plotSt.figPar(ff).axesOrd = 1:prod(plotSt.figPar(ff).axesNum);
    plotSt.figPar(ff).cols = cbrewer('qual', 'Paired', 2*max(plotSt.figPar(ff).relNumCol,3)); % so that it wont report a warning
end


%% managing the plotParSt input
if nargin == 5
    assert(isfield(plotParSt, 'figPar'), 'plotParSt is missing figPar field')
    
    % if input structure has length one, apply parameters on all figures
    if length(plotParSt.figPar) == 1 && numFigs > 1
        for ff=2:numFigs
            plotParSt.figPar(ff) = plotParSt.figPar(1);
        end
    end
    
    for ff=1:numFigs
        if isfield(plotParSt.figPar, 'axesNum')
            if ~isempty(plotParSt.figPar(ff).axesNum)
                tempAxeNum = plotParSt.figPar(ff).axesNum;
                assert(length(tempAxeNum) == 2, 'axesNum input for fig %d has wrong dimension', ff)
                assert(prod(tempAxeNum) >= plotSt.figPar(ff).relNumAxeRow * plotSt.figPar(ff).relNumAxeCol, 'axesNum for fig %d is too small', ff)
                plotSt.figPar(ff).axesNum = tempAxeNum;
            end
        end
        
        if isfield(plotParSt.figPar, 'axesOrd')
            if ~isempty(plotParSt.figPar(ff).axesOrd)
                tempAxeOrd = plotParSt.figPar(ff).axesOrd;
                assert(length(tempAxeOrd) == plotSt.figPar(ff).relNumAxeRow * plotSt.figPar(ff).relNumAxeCol, ... 
                       'axes order for fig %d is not equal to number of axes', ff) % assumes the matrix of rowXcol is full - maybe change
                plotSt.figPar(ff).axesOrd = tempAxeOrd;
            end
        end
        
        if isfield(plotParSt.figPar, 'cols')
            if ~isempty(plotParSt.figPar(ff).cols)
                tempCols = plotParSt.figPar(ff).cols;
                assert(size(tempCols, 2) == 3, 'Color input should be a 3 element vector')
                assert(size(tempCols, 1) >= plotSt.figPar(ff).relNumCol, 'Colors in fig %d exceed colors inputs', ff)
                plotSt.figPar(ff).cols = tempCols;
            end
        end
        
    end
end
   
%% plotting the data

maxTime = 0;
sbCol = [1,1,1]*0.6; 
fbCol = [1,1,1]*0.8; 
stimY = -7.5; 

for ff=1:numFigs
    
    handles.fig(ff).figH = figure;
    relAxe = plotSt.figPar(ff).axesNum;
    relOrd = plotSt.figPar(ff).axesOrd;
    relCol = plotSt.figPar(ff).cols;
    posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.01, 0.01, relAxe);
    tempFigInd = relTable{:, figName} == figIndVals(ff);
    
    axeCount = 0;
    for ii=1:plotSt.figPar(ff).relNumAxeRow
        
        tempAxeRowInd = relTable{:, axeNames{1}} == plotSt.figPar(ff).relAxeRowVal(ii);
        tempAxeRowVal = plotSt.figPar(ff).relAxeRowVal(ii);
        
        for jj=1:plotSt.figPar(ff).relNumAxeCol
            
            axeCount= axeCount+1;
            tempAxeColInd = relTable{:, axeNames{2}} == plotSt.figPar(ff).relAxeColVal(jj);
            tempAxeColVal = plotSt.figPar(ff).relAxeColVal(jj);
            
            handles.fig(ff).axh(ii,jj) = axes('position', posCell{relOrd(axeCount)});
            hold on
            
            linColInd = 1:2:2*plotSt.figPar(ff).relNumCol;
            datColInd = 2:2:2*plotSt.figPar(ff).relNumCol;
        
            for kk=1:plotSt.figPar(ff).relNumCol
                
                tempColInd = relTable{:, colName} == plotSt.figPar(ff).relColVal(kk);
            
                rowInd = prod([tempFigInd, tempAxeRowInd, tempAxeColInd, tempColInd], 2);
                gratInd = relTable{logical(rowInd), 'index'};
                
                if isempty(gratInd)
                    continue
                end
                
                assert(length(gratInd) == 1, 'variable names combination does not provide unique identity in gratingTable')
            
%                 plotIndsSt = getStimInds(mmResStruct, [gratInd, nan, nan, nan]);
%                 plotInds = [plotIndsSt(:).inds];
%                 [datToPlot, timeToPlot] = getStimDataByInds(mmResStruct, plotInds);
                
                datToPlot = mmResStruct.mmResult(gratInd).subData.baseSub(:,2);
                timeToPlot = mmResStruct.mmResult(gratInd).subData.baseSub(:,1);
                sbInd = mmResStruct.mmResult(gratInd).sbInd; 
                
                
                if ~isempty(mmResStruct.mmResult(gratInd).linSum)
                    linToPlot = mmResStruct.mmResult(gratInd).linSum(:,2);
                    plot(timeToPlot, linToPlot, 'lineWidth', 2, 'color', relCol(linColInd(kk), :))
                end
                
                % with the way color are processed here it will not work with
                % regular color input
                
                notNanInds = all(~isnan(timeToPlot), 2);
                meanTime = mean(timeToPlot(notNanInds, :),2);
                
                plot(timeToPlot, datToPlot, 'lineWidth', 2, 'color', relCol(datColInd(kk), :))
                titString = sprintf('%s: %d %s: %d', axeNames{1}, allAx1Vals(logical(rowInd)),axeNames{2}, allAx2Vals(logical(rowInd))); 
                title(titString, 'fontsize', 10)
                
                if ~isempty(sbInd)
                    line(timeToPlot(sbInd), [stimY, stimY], 'linewidth', 4, 'color', sbCol)
                    line([0, diff(timeToPlot(sbInd))], [stimY, stimY]-0.5, 'linewidth', 4, 'color', fbCol)
                end
                
                if ii==1
                    title(num2str(tempAxeColVal))
                end
                
                if jj==1
                    ylabel(num2str(tempAxeRowVal))
                end
                
                if meanTime(end) > maxTime
                    maxTime =  meanTime(end);
                end
            end
            
        end
        
    end
    
end



yyMax = cell(1,numFigs); 
yyMin = cell(1,numFigs); 

for ff=1:numFigs
    
    allLinesH = findobj(handles.fig(ff).axh, 'Type', 'Line', '-and', 'lineWidth', 2);
    
    for ii=1:length(allLinesH)
        
        yDat = allLinesH(ii).YData; 
        yyMax{ff} = [yyMax{ff} , max(yDat)];
        yyMin{ff} = [yyMin{ff} , min(yDat)];
        
    end
end

totMax = cellfun(@max, yyMax); 
totMin = cellfun(@min, yyMin); 


for ff=1:numFigs
    yBuf = (totMax(ff) - totMin(ff))/10; 
    totYLim = [totMin(ff) - yBuf, totMax(ff) + yBuf];
    for ii=1:size(handles.fig(ff).axh, 1)
        for jj=1:size(handles.fig(ff).axh, 2)
            handles.fig(ff).axh(ii,jj).YLim = totYLim;
            handles.fig(ff).axh(ii,jj).XLim = [-100, maxTime];
        end
    end
    
    midRow = ceil(ii/2);
    midCol = ceil(jj/2);
%     xxTick = handles.fig(ff).axh(1,1).XTick;
%     yyTick = handles.fig(ff).axh(1,1).YTick;
%     xxTickLab = handles.fig(ff).axh(1,1).XTickLabel;
%     yyTickLab = handles.fig(ff).axh(1,1).YTickLabel;
    
    set(handles.fig(ff).axh(:, 2:end), 'YTickLabel', [])
    set(handles.fig(ff).axh(1:end-1, :), 'XTickLabel', [])
    
    tempYLab = handles.fig(ff).axh(midRow, 1).YLabel.String;
    handles.fig(ff).axh(midRow, 1).YLabel.String = {axeNames{1}; tempYLab};
    
    tempXTit = handles.fig(ff).axh(1, midCol).Title.String;
    handles.fig(ff).axh(1, midCol).Title.String = {axeNames{2}; tempXTit};
    
end

        
if nargout == 1
    varargout{1} = handles;
end
    

end