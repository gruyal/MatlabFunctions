function varargout = plotStaticGratingDataandModel(allCellDat)

% function varargout = plotStaticGratingDataandModel(allCellDat)
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



sgTableInd = allCellDat.table.protType == 2;
sgTable = allCellDat.table(sgTableInd, :);

allPhase = unique(sgTable.phase);
relDur = [40, 160];
xxLim = [-50, 750];

posCell = generatePositionCell(0.05, 0.975, 0.075, 0.95, 0.01, 0.05, [length(allPhase), length(relDur)]);
axh = gobjects(size(posCell)); 

preCol = cbrewer('div', 'BrBG', 8); 
pCol = preCol([3,2,1], :);
modCol = [1,1,1]*0.8;


titText = {'stepDur 40ms', 'relDur 160ms'};


figure('position', [900, 900, 1250, 300])

for dd=1:length(relDur)

    for pp=1:length(allPhase)

        relInd = sgTable.index(sgTable.phase == allPhase(pp) & sgTable.duration == relDur(dd));  
        
        
        dataV = allCellDat.SG(relInd).data;
        modelV = allCellDat.SG(relInd).model;
        timeV = allCellDat.SG(relInd).time;

        axh(pp, dd) = axes('position', posCell{pp, dd}); 
        hold on 

        plot(timeV, modelV, 'linewidth', 2, 'color', modCol)
        plot(timeV, dataV, 'linewidth', 2, 'color', pCol(2, :))


        hold off

        axh(pp,dd).XLim = xxLim; 

        if pp==1
            axh(pp,dd).YLabel.String = titText{dd};
            axh(pp,dd).YLabel.FontWeight = 'bold';
        end

        if dd < length(relDur)
            axh(pp,dd).XColor = 'none';
        end

    end

end

for jj=1:length(allPhase)
    axh(jj,1).Title.String = ['Phase:', num2str(jj)]; 
end

for dd = 1:length(relDur)
    tmpY = vertcat(axh(:,dd).YLim);
    for jj = 1:length(allPhase)
        axh(jj,dd).YLim = [min(tmpY(:,1)), max(tmpY(:,2))];
        if jj > 1
            axh(jj,dd).YColor = 'none';
        end
    end 
end
    

if nargout == 1
    varargout{1} = axh;
end

    
end
            
            
            
            