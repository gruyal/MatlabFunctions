function varargout = plotResultsFromrevPhiRotV7AFixed(processedFileDir, saveFlag)

% function plotResultsFromrevPhiRotV7AFixed(processedFileDir)
%
% this function plots the results for the processed data in processedFileName
% it is specifically designed for the experiment in revPhiRotationV7A06-07-21_10-53-05
%
% after plotting if saveFlag is true also saves the flies in the given
% directory


if nargin < 2
    saveFlag =0;
end

load(fullfile(processedFileDir, '/processedData.mat'), 'timeseries', 'timestamps', 'LmR_avg_over_reps', 'interhistogram', 'channelNames');


%% fixed pattern V7A structure

uSD = [20, 40, 80, 160, 320]; 
uWid = [2,5];

index = (1:40)';
cwFlag = repmat([1;1;0;0], 10, 1);
revPhiFlag = repmat([0;1], 20, 1);
stepDur = reshape(repmat(uSD, 8, 1), [],1); 
barWid = repmat(reshape(repmat(uWid, 4, 1), [],1), 5, 1); 

patTab = table(index, barWid, cwFlag, revPhiFlag, stepDur);

%% plotting LmR results per trial 

xTit = {'CW W2 S'; ''; 'CCW W2 S'; ''; 'CW W5 S'; ''; 'CCW W5 S'; ''}; 
yTit = arrayfun(@num2str, uSD, 'UniformOutput', false);

numR = 5;
numC = 8;

lmrChI = find(strcmp(channelNames.timeseries, 'LmR'));
freqChI = find(strcmp(channelNames.timeseries, 'F_chan'));

relChI = [lmrChI, freqChI];

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.01, 0.01, [numC, numR]);
axh = gobjects(size(posCell));

numStim = height(patTab);
numReps = size(timeseries, 3);

pCol = cbrewer('qual', 'Paired', 6);
yyLim = [-5, 5; 0, 2.5];


for ch=1:length(relChI) % lmr and freq
    
    tempFH = figure('position', [220, 200, 1400, 750]);
    
    for ii=1:numStim

        relDat = (squeeze(timeseries(relChI(ch), ii, :, :)))';
        axh(ii) = axes('position', posCell{ii});
        hold on
        revF = patTab.revPhiFlag(ii);

        if revF == 0
            relC = pCol(1:2, :);
        else
            relC = pCol(5:6, :);
        end

        line([0, 2.2], [0,0], 'color', [1,1,1]*0.8)
        plot(timestamps, relDat, 'linewidth', 1, 'color', relC(1,:))
        if relChI == 1
            plot(timestamps, LmR_avg_over_reps(ii, :), 'linewidth', 2, 'color', relC(2,:))
        else
            plot(timestamps, mean(relDat, 2), 'linewidth', 2, 'color', relC(2,:))
        end

        if ii<=numC 
            title(xTit{ii})
        end

        if ~ismember(ii, 1:numC:numStim)
            axh(ii).YColor = 'none';
        else
            axh(ii).YLabel.String = yTit{floor(ii/numC)+1};
        end

        if ii<numStim-numC+1
            axh(ii).XColor = 'none';
        end

        hold off
        
        axh(ii).YLim = yyLim(ch, :);
    end
    
    if saveFlag
        print2PDF(fullfile(processedFileDir, ['resultsPerStimCh', num2str(relChI(ch))]), tempFH) 
    end
    
end

%% histogram for closed loop and LmR symmetry

tempFH2 = figure('position', [400, 600, 1000, 275]);

pCol = cbrewer('seq', 'YlOrRd', 9);

numPos = size(interhistogram, 2);
numTr = size(interhistogram, 1);
shiftP = round(numPos/2); 


prePosCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -0.2, 0.15, 2);
posCell = generatePositionCell(prePosCell{1}(1), sum(prePosCell{1}([1,3])), ...
                                prePosCell{1}(2), sum(prePosCell{1}([2,4])), 0.1, -0.01, 1+length(relChI));
axh = gobjects(size(posCell));

axh(1) = axes('position', posCell{1});
hold on 

tempHistD = circshift(sum(interhistogram(1:round(numTr/2), :)), shiftP);
[maxHDV1, maxHDP] = max(tempHistD);
hdFWHM = find(tempHistD > maxHDV1/2, 1, 'last') - find(tempHistD > maxHDV1/2, 1, 'first');
hdmaxPos = maxHDP - shiftP;

report.FWHM(1) = hdFWHM;
report.maxPos(1) = hdmaxPos;

plot(1:numPos, tempHistD, ...
     'linewidth', 2, 'color', pCol(3, :))
line([find(tempHistD > maxHDV1/2, 1, 'last'), find(tempHistD > maxHDV1/2, 1, 'first')], ...
     [maxHDV1/2, maxHDV1/2], 'color', pCol(3, :), 'linewidth', 4)

 
tempHistD = circshift(sum(interhistogram(round(numTr/2):end, :)), shiftP);
[maxHDV2, maxHDP] = max(tempHistD);
hdFWHM = find(tempHistD > maxHDV2/2, 1, 'last') - find(tempHistD > maxHDV2/2, 1, 'first');
hdmaxPos = maxHDP - shiftP;

report.FWHM(2) = hdFWHM;
report.maxPos(2) = hdmaxPos; 
 
plot(1:numPos, tempHistD, ...
    'linewidth', 2, 'color', pCol(5, :))
line([find(tempHistD > maxHDV2/2, 1, 'last'), find(tempHistD > maxHDV2/2, 1, 'first')], ...
     [maxHDV2/2, maxHDV2/2], 'color', pCol(5, :), 'linewidth', 4)
 
text(20, max(maxHDV1, maxHDV2) * 9/10, ['FWHM: ', num2str(report.FWHM)])
title('inter-trial fixation')
hold off

for ch=1:length(relChI)

    axh(ch+1) = axes('position', posCell{ch+1});
    hold on 
    for ii=1:numReps
        allTSLmR = squeeze(timeseries(relChI(ch), :, ii, :));
        allTSLmR = allTSLmR(:);
        [binCount, binEdge] = histcounts(allTSLmR, 100);
        binCen = mean([binEdge(1:end-1); binEdge(2:end)]);

        plot(binCen, binCount, 'Color', pCol(ii+2, :), 'linewidth', 2)
        title(['intra-trial value ch', num2str(relChI(ch))]) 

    end

    hold off
    
    if ch==length(relChI)
        legend(arrayfun(@(x) ['repeat', num2str(x)], 1:5, 'uniformoutput', 0), 'location', 'northwest')
    end


end

numBins = 100; 
LmRRange = [-5, 5];

posCell2 = generatePositionCell(prePosCell{2}(1), sum(prePosCell{2}([1,3])), ...
                                prePosCell{2}(2), sum(prePosCell{2}([2,4])), 0.05, -0.01, numReps);
axh3 = gobjects(size(posCell2));

for ii=1:numReps
    allTSLmR = squeeze(timeseries(lmrChI, :, ii, :));
    allTSLmR = allTSLmR(:);
    binEdge = linspace(LmRRange(1), LmRRange(2), numBins+1);
    [binCount, binEdge] = histcounts(allTSLmR, binEdge);
    binCen = mean([binEdge(1:end-1); binEdge(2:end)]);
    
    % folding the histogram
    cenPosI = ceil(numBins/2);

    lhCount = fliplr(binCount(1:cenPosI));
    lhCen = -fliplr(binCen(1:cenPosI));

    rhCount = binCount(cenPosI+1:end);
    rhCen = binCen(cenPosI+1:end);
    
    totCount = sum([rhCount, lhCount]);
    normRHCount = rhCount / totCount;
    normLHCount = lhCount / totCount;
    
    report.symInd(ii) = mean(abs((normRHCount - normLHCount)) ./ mean(normRHCount + normLHCount));
    
    axh3(ii) = axes('position', posCell2{ii});
    hold on 
    plot(lhCen, normLHCount, 'Color', pCol(ii+2, :), 'linewidth', 1)
    plot(rhCen, normRHCount, 'Color', pCol(ii+2, :), 'linewidth', 2)
    title(['LmR sym index', num2str(report.symInd(ii))]) 

end

hold off

if nargout > 0
    varargout{1} = report;
end


if saveFlag
    print2PDF(fullfile(processedFileDir, ['resultsHistograms', num2str(relChI(ch))]), tempFH2) 
end

    
end

