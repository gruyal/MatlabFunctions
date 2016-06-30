function varargout = plotMinMovLinSumResults(minMovLinSumStruct)

% function plotMinMovLinSumResults(minMovLimSumStruct)
%
% This function plots the results from calcMinMovingBarBasedOnSingleBar. 
% input strcture should be the output from the above function

datSiz = size(minMovLinSumStruct);

axh = zeros(datSiz);
posCell = generatePositionCell(0.05, 0.975, 0.025, 0.975, 0.02, 0.02, [4, datSiz(2)]);

pCol = cbrewer('qual', 'Paired', 4);

for ii=1:datSiz(1)
    
    stP = minMovLinSumStruct(ii,1,2).data.table.startPos;
    endP = minMovLinSumStruct(ii,1,2).data.table.stopPos;
    
    figure('name', ['positions', num2str(stP), ':', num2str(endP)])
    

    
    for jj=1:datSiz(2)
        
        axh(ii,jj,1) = axes('position', posCell{1, jj}); % comparing 2 directions
        
        plot(minMovLinSumStruct(ii,jj,1).subData.baseSub(:,1), minMovLinSumStruct(ii,jj,1).subData.baseSub(:,2), ...
             'linewidth', 3, 'color', pCol(2,:))
        hold on
        plot(minMovLinSumStruct(ii,jj,2).subData.baseSub(:,1), minMovLinSumStruct(ii,jj,2).subData.baseSub(:,2), ...
             'linewidth', 3, 'color', pCol(4,:))
         
        yLab = get(gca, 'ylabel');
        set(yLab, 'string', ['stepDur:', num2str(minMovLinSumStruct(ii,jj,1).data.table.stepDur)])
        
        xMin = min(minMovLinSumStruct(ii,jj,1).subData.baseSub(1,1), minMovLinSumStruct(ii,jj,2).subData.baseSub(1,1));
        xMax = max(minMovLinSumStruct(ii,jj,1).subData.baseSub(end,1), minMovLinSumStruct(ii,jj,2).subData.baseSub(end,1));
         
        axh(ii,jj,2) = axes('position', posCell{2, jj}); % comparing 1 dir to linSum
        
        plot(minMovLinSumStruct(ii,jj,1).subData.baseSub(:,1), minMovLinSumStruct(ii,jj,1).subData.baseSub(:,2), ...
             'linewidth', 3, 'color', pCol(2,:))
        hold on 
        plot(minMovLinSumStruct(ii,jj,1).linSum(:,1), minMovLinSumStruct(ii,jj,1).linSum(:,2), ...
             'linewidth', 3, 'color', pCol(1,:))
         
        axh(ii,jj,3) = axes('position', posCell{3, jj}); % comparing 2 dir to linSum
        
        plot(minMovLinSumStruct(ii,jj,2).subData.baseSub(:,1), minMovLinSumStruct(ii,jj,2).subData.baseSub(:,2), ...
             'linewidth', 3, 'color', pCol(4,:))
        hold on 
        plot(minMovLinSumStruct(ii,jj,2).linSum(:,1), minMovLinSumStruct(ii,jj,2).linSum(:,2), ...
             'linewidth', 3, 'color', pCol(3,:))
         
        axh(ii,jj,4) = axes('position', posCell{4, jj}); % comparing 2 dir to linSum
        
         plot(minMovLinSumStruct(ii,jj,1).linSum(:,1), minMovLinSumStruct(ii,jj,1).linSum(:,2), ...
             'linewidth', 3, 'color', pCol(1,:))
        hold on 
        plot(minMovLinSumStruct(ii,jj,2).linSum(:,1), minMovLinSumStruct(ii,jj,2).linSum(:,2), ...
             'linewidth', 3, 'color', pCol(3,:))
         
        set(axh(ii,jj,:), 'xlim', [xMin, xMax])
        yyLim = get(axh(ii,jj,:), 'ylim');
        yyLim = vertcat(yyLim{:});
        
        set(axh(ii,jj,:), 'ylim', [min(yyLim(:,1)), max(yyLim(:,2))])
         
    end
    
    
    
end

set(axh(:,:, 2:4), 'yticklabel', {})
set(axh(:,1:end-1, :), 'xticklabel', {})


if nargout ==1
    
    varargout{1} = axh;
    
end


end



