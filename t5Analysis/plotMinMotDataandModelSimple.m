function varargout = plotMinMotDataandModelSimple(allCellDat)

% function plotMinMotDataandModelSimple(allCellDat)
%
% A modification of plotMinMotDataandModel that does not look at the
% gratingTable for minMot (since it is absent from the new version of organizingClusterData
% and simply plot data vs model for all stimuli present
%
% compares modeling results for all minimal motion combinations. 
% PD is red  ND blue and single bar (diagonal) grey.
% first and second peaks are also labeled for both data (white) and
% model(black)
% grey bar underneath plot is stimulus timing (for first and second bar)
% width corresponds to bar width (only 2, 4 are shown for T5 cells)
%
% INPUT
% allCellDat -          generated using organizingClusterData for a specific cell and iteration 
%
% OUTPUT 
% axh -                 optional. axes handles



close all

figure('position', [800, 20, 1500, 1300])

numPlotsPerR = 10; 
relDur = [40, 160];
relW = [2, 4];

preTab = allCellDat.table; 
mmTable = preTab(preTab.protType == 4, :); 


relTab = mmTable(ismember(mmTable.width, relW)  ... 
                 & ismember(mmTable.duration, relDur) ...
                 , :); 

numStim = height(relTab);
numC = ceil(numStim/numPlotsPerR); 

posCell = generatePositionCell(0.025, 0.975, 0.025, 0.975, 0.005, 0.005, [numPlotsPerR, numC]);
axh = gobjects(size(posCell)); 
preCol = cbrewer('qual', 'Paired', 6);
pColD = preCol([2,6], :);
pColD(3,:) = [1,1,1]*0.65;
pColM = preCol([1,5], :);
pColM(3,:) = [1,1,1]*0.8;
fontW = {'normal'; 'bold'};
stimY = -2; 

for ii=1:numStim

    relInd = relTab.index(ii); 
    tempSt = allCellDat.MM(relInd);
    tempDat = tempSt.data;
    tempMod = tempSt.model; 
    tempT = tempSt.time; 
    tempDir = relTab.direction(ii); 

    datMax = relTab.maxData(ii,:); 
    modMax = relTab.maxMod(ii,:); 
    durI = relTab.duration(ii);
    wid = relTab.width(ii);


    if tempDir == 1
        colInd = 2;
    elseif tempDir == 0
        colInd = 1;
    else
        colInd = 3; % for the same position
    end

    axh(ii) = axes('position', posCell{ii}); 
    hold on 
    
    plot(tempT, tempDat, 'linewidth', wid, 'color', pColD(colInd, :))
    plot(tempT, tempMod, 'linewidth', wid, 'color', pColM(colInd, :))
    
    line([0, durI], [-2, -2], 'linewidth', 3, 'color', pColM(3, :))

    
    hold off

end

axh = axh(1:numStim);

allYLim = vertcat(axh(:).YLim);
globYLim = [min(allYLim(:,1)), max(allYLim(:,2))];

for ii=1:numel(axh)
    axh(ii).YLim = globYLim;
    axh(ii).YColor = 'none';
    axh(ii).XColor = 'none';

end

if nargout == 1
    varargout{1} = axh;
end


    
end

