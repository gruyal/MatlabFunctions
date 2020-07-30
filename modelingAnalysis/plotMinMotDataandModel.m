function varargout = plotMinMotDataandModel(allCellDat)

% function plotMinMotDataandModel(allCellDat)
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
relWandSP = [2,1; 4,0];

preTab = allCellDat.table; 
mmTable = preTab(preTab.protType == 4, :); 


relTab = mmTable((mmTable.width == relWandSP(1,1) & mmTable.speedCor == relWandSP(1,2) | ...
                 mmTable.width == relWandSP(2,1) & mmTable.speedCor == relWandSP(2,2)) ...
                 & ismember(mmTable.duration, relDur) ...
                 & mmTable.timeDiff > 0 ...
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

    fbp = relTab.FBPos(ii);
    sbp = relTab.SBPos(ii);
    fbTime = relTab.fbTime(ii, :);
    sbTime = relTab.sbTime(ii, :);
    datMax = relTab.maxData(ii,:); 
    modMax = relTab.maxMod(ii,:); 
    durI = double(relTab.duration(ii) == 160) + 1;
    wid = relTab.width(ii);


    if sbp > fbp
        colInd = 2;
    elseif sbp < fbp
        colInd = 1;
    else
        colInd = 3; % for the same position
    end

    axh(ii) = axes('position', posCell{ii}); 
    hold on 

    line([fbTime, nan, sbTime], [1,1,nan,1.5,1.5]*stimY, 'linewidth', wid, 'color', pColD(3,:)); 
    plot(tempT, tempDat, 'linewidth', 2, 'color', pColD(colInd, :))
    plot(tempT, tempMod, 'linewidth', 2, 'color', pColM(colInd, :))
    plot([mean(fbTime), mean(sbTime)], datMax, '-o', 'linewidth', wid/2, 'color', pColD(colInd, :), 'markerfacecolor', 'w', 'markeredgecolor', 'k')
    plot([mean(fbTime), mean(sbTime)], modMax, '-o', 'linewidth', wid/2, 'color', pColM(colInd, :), 'markerfacecolor', 'k', 'markeredgecolor', 'none')
    title([num2str(fbp), 'to',  num2str(sbp)], 'FontWeight', fontW{durI})
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

