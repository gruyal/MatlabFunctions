
function varargout = plotbaseSubSingleBar(singleBarSt)

% plotbaseSubSingleBar(singleBarSt)
%
% this function takes the output of generateAlignedSingleBarSt and plots it
% together with response estimate



xRange = [-200, 600; ...
          -200, 600; ... 
          -200, 800; ...
          -200, 800];

datSiz = size(singleBarSt);

posCell = generatePositionCell(0.05, 0.975, 0.025, 0.975, 0.005, 0.005, [datSiz(1), datSiz(2)]);

axh = zeros(datSiz);

pCol = cbrewer('qual', 'Paired', 2*datSiz(2));

figure

for ii=1:datSiz(1)
    
    for jj=1:datSiz(2)
        
        axh(ii,jj) = axes('position', posCell{ii,jj});
        
        plot(singleBarSt(ii,jj).subData.baseSub(:,1), singleBarSt(ii,jj).subData.baseSub(:,2), 'linewidth', 3, 'color', pCol(2*jj,:))
        hold on 
        plot(singleBarSt(ii,jj).subData.baseSub(:,1), singleBarSt(ii,jj).subData.baseSubMed, 'linewidth', 3, 'color', pCol(2*jj-1,:))
        
        line(xRange(jj, :), [singleBarSt(ii,jj).maxResp, singleBarSt(ii,jj).maxResp], 'color', pCol(2*jj-1, :), 'linewidth', 3)
        line(xRange(jj, :), [singleBarSt(ii,jj).minResp, singleBarSt(ii,jj).minResp], 'color', pCol(2*jj-1, :), 'linewidth', 2)
        line(xRange(jj, :), [0, 0], 'color', [1,1,1]*0.8, 'linewidth', 2)
        hold off
    end
    
end

for ii=1:datSiz(2) 
    yyLim = get(axh(:,ii), 'ylim');
    yyLim = vertcat(yyLim{:});
    set(axh(:,ii), 'ylim', [min(yyLim(:,1)), max(yyLim(:,2))])
    set(axh(:,ii), 'xlim', xRange(ii,:))
end

set(axh(2:end, :), 'yticklabel', {})
set(axh(:, 1:end-1), 'xticklabel', {})


if nargout==1
    varargout{1} = axh;
end





end