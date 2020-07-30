function varargout = plotMovingBarandModel(allCellDat)

% function varargout = plotMovingBarandModel(allCellDat)
%
% This function plots data and model comparisons for PD and ND in relevant
% durations only
%
% INPUT 
% allCellDat -      generated using organizingClusterData for a specific
%                   cell and iteration 
%
% OUTPUT
% axh               optional. axes handles 




relWid = [1,2,4];
relDur = [40, 160];
xxLim = [-50, 1000; -50, 2200];

posCell = generatePositionCell(0.075, 0.95, 0.05, 0.95, 0.05, 0.05, [length(relDur), length(relWid)]);
axh = gobjects(size(posCell)); 

preCol = cbrewer('qual', 'Paired', 6); 
pCol = preCol([1,2,5,6], :);
titText = {'stepDur 40ms', 'relDur 160ms'};

mbInds = allCellDat.table.protType == 3;

mbTable = allCellDat.table(mbInds, :); 

figure('position', [1600, 440, 650, 750])

for ww=1:length(relWid)

    for dd=1:length(relDur)

        relInds = mbTable.index(mbTable.width == relWid(ww) & mbTable.duration == relDur(dd));  
        relDirs = mbTable.direction( mbTable.width == relWid(ww) & mbTable.duration == relDur(dd));
  
        if isempty(relInds)
            continue
        end
        
        assert(length(relDirs) == 2, 'missing direction')
        assert(relDirs(1) == 0 & relDirs(2) == 1, 'directions are flipped')
        
        dataND = allCellDat.MB(relInds(1)).data;
        modelND = allCellDat.MB(relInds(1)).model;
        timeND = allCellDat.MB(relInds(1)).time;

        dataPD = allCellDat.MB(relInds(2)).data;
        modelPD = allCellDat.MB(relInds(2)).model;
        timePD = allCellDat.MB(relInds(2)).time;

        axh(dd, ww) = axes('position', posCell{dd, ww}); 
        hold on 

        plot(timeND, modelND, 'linewidth', 2, 'color', pCol(1, :))
        plot(timeND, dataND, 'linewidth', 2, 'color', pCol(2, :))

        plot(timePD, modelPD, 'linewidth', 2, 'color', pCol(3, :))
        plot(timePD, dataPD, 'linewidth', 2, 'color', pCol(4, :))

        hold off

        axh(dd,ww).XLim = xxLim(dd,:); 

        if dd==1
            axh(dd,ww).YLabel.String = ['Width:', num2str(relWid(ww))];
        end

        if ww < length(relWid)
            axh(dd,ww).XColor = 'none';
        end

    end

end

for jj=1:length(relDur)
    axh(jj,1).Title.String = titText(jj); 
end

axh = findobj(axh, 'xcolor', 'none'); % in case some axes are empty

legend(axh(1, end), 'ModND', 'DatND', 'ModPD', 'DatPD')
legend(axh(1, end), 'boxoff')

if nargout == 1
    varargout{1} = axh;
end

    
end
            
            
            
            