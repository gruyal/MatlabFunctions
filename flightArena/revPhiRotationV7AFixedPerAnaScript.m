

% relDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/flightArena/revPhiExp/revPhiRotationV7A06-07-21_10-53-05';
relDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/flightArena/revPhiExp/revPhiRotationV7B06-07-21_12-54-02';


load(fullfile(relDir, '06_14_2021/SS2344 X JFRC28-11_07_27/processedData.mat'));

%% fixed pattern V7A structure

uSD = [20, 40, 80, 160, 320]; 
uWid = [2,5];

index = (1:40)';
cwFlag = repmat([1;1;0;0], 10, 1);
revPhiFlag = repmat([0;1], 20, 1);
stepDur = reshape(repmat(uSD, 8, 1), [],1); 
barWid = repmat(reshape(repmat(uWid, 4, 1), [],1), 5, 1); 

patTab = table(index, barWid, cwFlag, revPhiFlag, stepDur);

clear index cwFlag stepDur revPhiFlag barWid

%%

xTit = {'CW W2 S'; ''; 'CCW W2 S'; ''; 'CW W5 S'; ''; 'CCW W5 S'; ''}; 
yTit = arrayfun(@num2str, uSD, 'UniformOutput', false);

numR = 5;
numC = 8;

lmrChI = 4;
posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.01, 0.01, [numC, numR]);
axh = gobjects(size(posCell));

numStim = height(patTab);
numReps = size(timeseries, 3);

pCol = cbrewer('qual', 'Paired', 6);
yyLim = [-5, 5];


for ii=1:numStim
    
    relDat = (squeeze(timeseries(lmrChI, ii, :, :)))';
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
    plot(timestamps, LmR_avg_over_reps(ii, :), 'linewidth', 2, 'color', relC(2,:))
    
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
    axh(ii).YLim = yyLim;
end
    
    
%% histogram for closed loop and LmR symmetry

pCol = cbrewer('seq', 'YlOrRd', 9);

numPos = size(interhistogram, 2);
numTr = size(interhistogram, 1);
shiftP = round(numPos/2); 

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.1, -0.01, 2);
axh = gobjects(size(posCell));

axh(1) = axes('position', posCell{1});
hold on 
plot(1:numPos, circshift(sum(interhistogram(1:round(numTr/2), :)), shiftP), ...
     'linewidth', 2, 'color', pCol(3, :))
plot(1:numPos, circshift(sum(interhistogram(round(numTr/2):end, :)), shiftP), ...
    'linewidth', 2, 'color', pCol(5, :))
hold off


axh(1) = axes('position', posCell{2});
hold on 
for ii=1:numReps
    allTSLmR = squeeze(timeseries(lmrChI, :, ii, :));
    allTSLmR = allTSLmR(:);
    [binCount, binEdge] = histcounts(allTSLmR, 100);
    binCen = mean([binEdge(1:end-1); binEdge(2:end)]);
        
    plot(binCen, binCount, 'Color', pCol(ii+2, :), 'linewidth', 2)
    
end

hold off


    
    
    