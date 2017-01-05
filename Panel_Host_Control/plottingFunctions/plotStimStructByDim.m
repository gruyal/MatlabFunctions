function varargout = plotStimStructByDim(pStruct, figDim, axeDim, colDim, plotParSt)

% function varargout = plotStimStructByDim(pStruct, figDim, posDim, colDim, plotSt)
%
% This functions plots the data in pStruct based on the separation in the
% different dimensions. pStruct should have a .stim field, and should work
% with getStimInds (have either a relInds or a combInds field). 
% 
% INPUTS
% pStruct -         protocol structure with .stim field and .data in it
% figDim -          dimension accorind to which data would be seperated
%                   into different figures (relInds or combInds).
% axeDim -          same as above only to different axes within figure.
% colDim -          same as above only to different colors within axes. 
%
% plotParSt -          (optional) structure with different parameters that
%                   will affect data plotting 
% .figPar -        fields within .figPar will be applied on each figure.
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
%   .plotReps -     logical flag for plotting repeats
%
% if figPar is left empty for a particular figure it will be generated by
% default parameters.


%% verifying input correctness

% close all


assert(isfield(pStruct, 'stim'), 'pStrcut is missing stim field')
assert(isfield(pStruct.stim, 'data'), 'stim is missing data field')

if isfield(pStruct.stim, 'relInds')
    allInds = vertcat(pStruct.stim.relInds);
elseif isfield(pStruct.stim, 'combInds')
    allInds = vertcat(pStruct.stim.combInds);
else
    error('stim field is missing relInds or combInds field')
end

indsLen = size(allInds,2);

assert(prod([figDim, axeDim, colDim] <= size(allInds, 2)) == 1, 'specified dimension exceed the length of stim indices')

figIndVals = unique(allInds(:, figDim));
numFigs = length(figIndVals);


for ii=1:numFigs
    axeIndVals = unique(allInds(allInds(:, figDim) == figIndVals(ii), axeDim));
    colIndVals = unique(allInds(allInds(:, figDim) == figIndVals(ii), colDim));
    plotSt.figPar(ii).relAxeVal = axeIndVals;
    plotSt.figPar(ii).relColVal = colIndVals;
    plotSt.figPar(ii).relNumAxe = length(axeIndVals);
    plotSt.figPar(ii).relNumCol = length(colIndVals);
    plotSt.figPar(ii).axesNum = [ceil(sqrt(plotSt.figPar(ii).relNumAxe)), ceil(sqrt(plotSt.figPar(ii).relNumAxe))];
    plotSt.figPar(ii).axesOrd = 1:length(axeIndVals);
    plotSt.figPar(ii).cols = cbrewer('qual', 'Paired', 2*plotSt.figPar(ii).relNumCol); 
    plotSt.figPar(ii).plotReps = 1;
end


%% managing the plotParSt input
if nargin == 5
    assert(isfield(plotParSt, 'figPar'), 'plotParSt is missing figPar field')
    
    % if input structure has length one, apply parameters on all figures
    if length(plotParSt.figPar) == 1 && numFigs > 1
        for ii=2:numFigs
            plotParSt.figPar(ii) = plotParSt.figPar(1);
        end
    end
    
    for ii=1:numFigs
        if isfield(plotParSt.figPar, 'axesNum')
            if ~isempty(plotParSt.figPar(ii).axesNum)
                tempAxeNum = plotParSt.figPar(ii).axesNum;
                assert(length(tempAxeNum) == 2, 'axesNum input for fig %d has wrong dimension', ii)
                assert(prod(tempAxeNum) >= plotSt.figPar(ii).relNumAxe, 'axesNum for fig %d is too small', ii)
                plotSt.figPar(ii).axesNum = tempAxeNum;
            end
        end
        
        if isfield(plotParSt.figPar, 'axesOrd')
            if ~isempty(plotParSt.figPar(ii).axesOrd)
                tempAxeOrd = plotParSt.figPar(ii).axesOrd;
                assert(length(tempAxeOrd) == plotSt.figPar(ii).relNumAxe, ...
                       'axes order for fig %d os not equal to number of axes', ii)
                plotSt.figPar(ii).axesOrd = tempAxeOrd;
            end
        end
        
        if isfield(plotParSt.figPar, 'cols')
            if ~isempty(plotParSt.figPar(ii).cols)
                tempCols = plotParSt.figPar(ii).cols;
                assert(size(tempCols, 2) == 3, 'Color input should be a 3 element vector')
                assert(size(tempCols, 1) >= plotSt.figPar(ii).relNumCol, 'Colors in fig %d exceed colors inputs', ii)
                plotSt.figPar(ii).cols = tempCols;
            end
        end
        
        if isfield(plotParSt.figPar, 'plotReps')
            if ~isempty(plotParSt.figPar(ii).plotReps)
                tempRepsFlag = plotParSt.figPar(ii).plotReps;
                assert(ismember(tempRepsFlag, [0,1]), 'plotReps should be logical ')
                plotSt.figPar(ii).plotReps = tempRepsFlag;
            end
        end
    end
end
   
%% plotting the data

maxTime = 0;

for ii=1:numFigs
    
    handles.fig(ii).figH = figure;
    relAxe = plotSt.figPar(ii).axesNum;
    relOrd = plotSt.figPar(ii).axesOrd;
    relCol = plotSt.figPar(ii).cols;
    posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.02, 0.02, relAxe);
    
    for jj=1:plotSt.figPar(ii).relNumAxe
        
        handles.fig(ii).axh(jj) = axes('position', posCell{relOrd(jj)});
        hold on
        
        if 2*plotSt.figPar(ii).relNumCol == size(plotSt.figPar(ii).cols, 1)
            indvColInd = 1:2:2*plotSt.figPar(ii).relNumCol;
            meanColInd = 2:2:2*plotSt.figPar(ii).relNumCol;
        else
            indvColInd = 1:plotSt.figPar(ii).relNumCol;
            meanColInd = 1:plotSt.figPar(ii).relNumCol;
        end
        
        for kk=1:plotSt.figPar(ii).relNumCol
            indsTemp = nan(1, indsLen);
            indsTemp([figDim, axeDim, colDim]) = [figIndVals(ii), plotSt.figPar(ii).relAxeVal(jj), plotSt.figPar(ii).relColVal(kk)];
            plotIndsSt = getStimInds(pStruct, indsTemp);
            plotInds = [plotIndsSt(:).inds];
            [datToPlot, timeToPlot] = getStimDataByInds(pStruct, plotInds);
            
            % with the way color are processed here it will not work with
            % regular color input
            if plotSt.figPar(ii).plotReps
                plot(timeToPlot, datToPlot, 'lineWidth', 1, 'color', relCol(indvColInd(kk), :))
            end
            
            notNanInds = all(~isnan(timeToPlot), 2);
            meanTime = mean(timeToPlot(notNanInds, :),2);
            
            plot(meanTime, mean(datToPlot(notNanInds, :),2), 'lineWidth', 3, 'color', relCol(meanColInd(kk), :))
            
            title(num2str(jj))
            if meanTime(end) > maxTime
                maxTime =  meanTime(end);
            end
            
        end
        
    end
    
end


allY = [];
for ii=1:numFigs
    allY = vertcat(allY, handles.fig(ii).axh(:).YLim);
end

totYLim = [min(allY(:,1)), max(allY(:,2))];

for ii=1:numFigs
    for jj=1:length(handles.fig(ii).axh)
        handles.fig(ii).axh(jj).YLim = totYLim;
        handles.fig(ii).axh(jj).XLim = [0, maxTime];
    end
end
        
if nargout == 1
    varargout{1} = handles;
end
    

end
