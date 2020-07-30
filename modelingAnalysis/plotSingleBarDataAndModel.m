function varargout = plotSingleBarDataAndModel(allCellDat)

% function varargout = plotSingleBarDataAndModel(allCellDat)
%
% This function plots the SPFR data overlaid by the model and the
% difference between Data and Model
% Note: for sake of clarity dies not plot all positions or all duration 
%
% INPUT
% allCellDat -      generated using orginizingClusterData function for a specific cell and iteration    



relWid = [1,2,4];
relDur = [40, 160];
relPos = -6:6; 
xxLim = [-50, 600];

posCell = generatePositionCell(0.025, 0.975, 0.025, 0.975, 0.001, 0.01, [length(relPos), length(relDur) * length(relWid)]);
axh = gobjects(size(posCell)); 

preCol = cbrewer('div', 'BrBG', 8); 
pCol = preCol([3,2,1], :);
modCol = [1,1,1]*0.8;

sbInds = allCellDat.table.protType == 0 ;

sbTable = allCellDat.table(sbInds,:); 

figure('position', [375, 250, 2000, 1000])
titFlag = false(length(relPos),1); 

for ww=1:length(relWid)

    for dd=1:length(relDur)

        yyLims = zeros(length(relPos), 2); 
        relPFlag = false(length(relPos),1); 

        for pp=1:length(relPos)

            relInd = sbTable.index(sbTable.width == relWid(ww) & sbTable.duration == relDur(dd) & sbTable.position == relPos(pp));  

            if isempty(relInd)
                continue
            end

            assert(length(relInd) == 1, 'index wrong')

            dataV = allCellDat.SB(relInd).data;
            modelV = allCellDat.SB(relInd).model;
            timeV = allCellDat.SB(relInd).time;

            axh(pp, 2*(ww-1)+dd) = axes('position', posCell{pp, 2*(ww-1)+dd}); 
            hold on 

            line(xxLim, [0, 0], 'linewidth', 1, 'color', modCol)
            plot(timeV, modelV, 'linewidth', 2, 'color', modCol)
            plot(timeV, dataV, 'linewidth', 2, 'color', pCol(ww, :))
            plot(timeV, dataV-modelV, 'linewidth', 1, 'color', 'k')

            hold off

            axh(pp, 2*(ww-1)+dd).XLim = xxLim; 
            axh(pp, 2*(ww-1)+dd).XColor = 'none';  


            yyLims(pp, :) = axh(pp, 2*(ww-1)+dd).YLim; 
            relPFlag(pp) = 1; 

            if ~titFlag(pp)
                title(num2str(relPos(pp)), 'fontsize', 14)
                titFlag(pp) = 1;
            end

        end

        % after I have all yyLims for this combination
        relYLims = yyLims(relPFlag, :); 

        yMin = min(relYLims(:,1)); 
        yMax = max(relYLims(:,2)); 

        rel2Pos = find(relPFlag); 

        for pp=1:length(rel2Pos)
            axh(rel2Pos(pp), 2*(ww-1)+dd).YLim = [yMin, yMax];
            if pp > 1
                axh(rel2Pos(pp), 2*(ww-1)+dd).YColor = 'none'; 
            else
                axh(rel2Pos(pp), 2*(ww-1)+dd).YLabel.String = ['Wid:', num2str(relWid(ww)), ' Dur:', num2str(relDur(dd))];
                axh(rel2Pos(pp), 2*(ww-1)+dd).YLabel.FontSize = 14;
            end


        end

    end

end

lH = findobj(gca, 'color', 'k');

legend(lH, 'Data - Model')
legend('boxoff')


if nargout == 1
    varargout{1} = axh;
end

    
    
end
            
            
            
            