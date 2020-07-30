
function varargout = plotBaseSubSingleBarDiffWandV(singleBarSt)

% plotbaseSubSingleBar(singleBarSt)
%
% this function takes the output of generateAlignedSingleBarStwMinMaxDiffWandV and plots it
% together with response estimate

close all
allStimDur = [0.02, 0.04, 0.08, 0.16, 0.32];

xRange = [-200, 1000; ...
          -200, 1000; ... 
          -200, 1200; ...
          -200, 1200; ...
          -200, 1200];

datSiz = size(singleBarSt); % dimensions are position, duration, width, and value
% to deal with singelton in the end (if there is only one value D or B)  
if length(datSiz) == 3
    datSiz(4) = 1; 
end

if isfield(singleBarSt, 'wPosMax')
    datSiz = datSiz - [1,1,0,0]; % since the last position for the first 2 dimensions is statistics
end



plotDur = zeros(1, datSiz(2));
plotPos = zeros(1, datSiz(1));

for jj=1:datSiz(2) %over all durations
    for ii=1:datSiz(1)
        for ww=1:datSiz(3)
            for vv=1:datSiz(4)
                if ~singleBarSt(ii,jj,ww,vv).empty
                    plotDur(jj) = singleBarSt(ii,jj,ww,vv).data.table.stimDur;
                    plotPos(ii) = singleBarSt(ii,jj,ww,vv).data.table.position;
                end
            end
        end
    end
end

plotDurInd = find(ismember(allStimDur, plotDur)); 


% axh = gobjects(datSiz(3), datSiz(1)-1, datSiz(2)-1); % -1 when you add
% the mormalized max and min 
axh = gobjects(datSiz(3), datSiz(1), datSiz(2));

for ww=1:datSiz(3) 
    
    posCell = generatePositionCell(0.05, 0.975, 0.025, 0.975, 0.005, 0.005, [datSiz(1), datSiz(2)]); % -1 
    
    axh(ww, :, :) = zeros(datSiz(1), datSiz(2)); %-1
    
    pCol = cbrewer('qual', 'Set1', 2*(datSiz(4)));
    markCol = cbrewer('qual', 'Paired', 2*(datSiz(4)));
    
    figure

    for ii=1:datSiz(1) % datSiz-1, since summary stat are in the empty end of the structure

        for jj=1:datSiz(2)

            axh(ww,ii,jj) = axes('position', posCell{ii,jj});
            hold on 
            
            for vv=1:datSiz(4)
                
                if singleBarSt(ii,jj,ww,vv).empty
                    continue
                end

                plot(singleBarSt(ii,jj,ww,vv).subData.baseSub(:,1), singleBarSt(ii,jj,ww,vv).subData.baseSub(:,2), 'linewidth', 4, 'color', pCol(vv,:))

                plot(singleBarSt(ii,jj,ww,vv).resp.maxTime, singleBarSt(ii,jj,ww,vv).resp.maxVal, 'o', ...
                    'markeredgecolor', 'r', 'markerfacecolor', markCol(2*vv, :), 'markersize', 8)

                plot(singleBarSt(ii,jj,ww,vv).resp.minTime, singleBarSt(ii,jj,ww,vv).resp.minVal, 'o', ...
                    'markeredgecolor', 'b', 'markerfacecolor', markCol(2*vv-1, :), 'markersize', 8)

                line(xRange(jj, :), [0, 0], 'color', [1,1,1]*0.8, 'linewidth', 2)
                
                
            end
            
            if jj==1
                title(num2str(plotPos(ii)))
            end
            
            hold off

            
        end

    end
        
end


ePosInd = singleBarSt(end, end, end, end).maxExtPosInd; 
iPosInd = singleBarSt(end, end, end, end).minInhPosInd; 


for ww=1:datSiz(3)
    for ii=1:datSiz(2)
        yyLim = vertcat(axh(ww, :,ii).YLim);
        for jj=1:datSiz(1)
            axh(ww,jj,ii).YLim = [min(yyLim(:,1)), max(yyLim(:,2))];
            axh(ww,jj,ii).XLim = xRange(plotDurInd(ii),:);
            
            if jj == ePosInd
                axh(ww,jj,ii).Box = 'on';
                axh(ww,jj,ii).YColor = 'r';
                axh(ww,jj,ii).XColor = 'r';
                axh(ww,jj,ii).YAxis.LineWidth = 2; 
                axh(ww,jj,ii).XAxis.LineWidth = 2; 
            end
            
            if jj == iPosInd
                axh(ww,jj,ii).Box = 'on';
                axh(ww,jj,ii).YColor = 'b';
                axh(ww,jj,ii).XColor = 'b';
                axh(ww,jj,ii).YAxis.LineWidth = 2; 
                axh(ww,jj,ii).XAxis.LineWidth = 2; 
            end
            
        end
    end
    
    set(axh(ww,2:end, :), 'yticklabel', {})
    set(axh(ww,:, 1:end-1), 'xticklabel', {})
end




if nargout==1
    varargout{1} = axh;
end





end