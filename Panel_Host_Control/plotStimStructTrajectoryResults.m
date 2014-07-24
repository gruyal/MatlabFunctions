function plotStimStructTrajectoryResults(pStruct, sepColDim)

% function plotStimStructResults(pStruct, sepColDim)
%
% This function reads the data from protocolStruct.stim and presents it
% separated by different figures (one for each direction) and colors based on the dimensions given
% (dimensions refers to the relInds component of the stim structure)

relChannel = 3; %since the first column in the data is timeStamp
% Choosing colors that will make sense (be able to see progression)
relCols = repmat(linspace(1, 0, 64), 3, 1)'; % gray scale
offset = 5;

% determining how to divide each figure (by checking maskPositions)
mPos = pStruct.maskPositions;
assert(iscell(mPos), 'mask positions should be a cell array')

newmPos = cellfun(@(x) sign(x(end,:)-x(1, :)), pStruct.maskPositions, 'uniformoutput', 0);
newmPos = vertcat(newmPos{:});
posTit = cellfun(@(x) [num2str(x(1, :)), ' to ', num2str(x(end, :))], pStruct.maskPositions, 'uniformoutput', 0);

allDirs = unique(newmPos, 'rows');
numDirs = size(allDirs, 1);
numTrajPerD = size(newmPos, 1)/numDirs;

dirInds = zeros(1, size(newmPos,1));
for ii=1:size(newmPos,1)
    tempLog = all(bsxfun(@eq, newmPos(ii,:), allDirs),2);
    dirInds(ii) = find(tempLog);
end

% xVals = unique(mPos(:,1));
% numX = length(xVals);
% yVals = unique(mPos(:,2));
% numY = length(yVals);

%flipping lr puts the arena lower position low in the subplots
subPlotPos = fliplr(generatePositionCell(0.05, 0.95, 0.05, 0.95, -1, 0.05, numTrajPerD));
allInds = vertcat(pStruct.stim.relInds);

% since the 4th dim is mask position which determines subplots
assert(ismember(sepColDim, 1:3), 'Dimension by which to designate colors cannot be bigger than 3')

numFigs = numDirs;
numCols = length(unique(allInds(:, sepColDim)));
goodColsInds = round(linspace(1+offset, size(relCols, 1) - offset, numCols)); 
axh = zeros(numFigs, numTrajPerD);


for ii=1:numFigs
    figure
    relDirInd = find(dirInds == ii);
    firstInds = cell(1, length(relDirInd));
    for ind=1:length(relDirInd)
        firstInds{ind} = allInds(:, 4) == relDirInd(ind);
    end
    set(gcf, 'name', ['Direction', num2str(allDirs(ii, :))])
    minX = 1000;
    maxX = 0;
    handForLegend = zeros(1, numCols);
    for jj=1:numTrajPerD
        axh(ii, jj) = axes('position', subPlotPos{jj});
        hold on 
        title(posTit{relDirInd(jj)})
        for kk=1:numCols
            secInds = allInds(:,sepColDim) == kk;
            plotInds = find(firstInds{jj}+secInds == 2);
            plotCol = relCols(goodColsInds(kk), :);
            for mm=1:length(plotInds)
                dataX = pStruct.stim(plotInds(mm)).data{1}(1, :); 
                if dataX(1) < minX
                    minX = dataX(1);
                end
                if dataX(end) > maxX
                    maxX = dataX(end);
                end
                dataY = pStruct.stim(plotInds(mm)).data{1}(relChannel, :)*10; % to convert to mV
                lh = plot(axh(ii, jj), dataX, dataY, 'linewidth', 1, 'color', plotCol);
            end
            handForLegend(kk) = lh;
        end
        hold off
    end
    allY = get(axh(ii, :), 'ylim');
    allY = [allY{:}];
    maxY = max(allY);
    minY = min(allY);
    set(axh(ii, :), 'ylim', [minY, maxY])
    set(axh(ii, :), 'xlim', [minX, maxX])
%     xlab = get(axh(ii, 1), 'xticklabel');
%     ylab = get(axh(ii, 1), 'yticklabel');
%     set(axh(ii, :), 'yticklabel', {}, 'xticklabel', {})
%     set(axh(ii, 1:numY:(numX*numY)), 'xticklabel', xlab)
%     set(axh(ii, 1:numY), 'yticklabel', ylab)
    
    legend(axh(ii, end), handForLegend, arrayfun(@num2str, 1:numCols, 'uniformoutput', 0))
end






end
