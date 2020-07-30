
function varargout = plotBaseSubSingleBarDiffWandVDiffMeans(singleBarSt1, singleBarSt2)

% plotbaseSubSingleBar(singleBarSt)
% 
% This function is a modification of plotBaseSubSingleBarDiffWandVDiff but
% is used to compare the 2 version of
% generateAlignedSingleBarStwMinMaxDiffWandV and thier effect on the mean
% 
%
% INPUTs
% singleBarSt1/2        - outputs from
%                         generateAlignedSingleBarStwMinMaxDiffWandV and
%                         generateAlignedSingleBarStwMinMaxDiffWandV2 

close all
allStimDur = [0.02, 0.04, 0.08, 0.16, 0.32];

xRange = [-200, 1000; ...
          -200, 1000; ... 
          -200, 1200; ...
          -200, 1200; ...
          -200, 1200];

datSiz = size(singleBarSt1); % dimensions are position, duration, width, and value
% to deal with singelton in the end (if there is only one value D or B)  
if length(datSiz) == 3
    datSiz(4) = 1; 
end

if isfield(singleBarSt1, 'wPosMax')
    datSiz = datSiz - [1,1,0,0]; % since the last position for the first 2 dimensions is statistics
end



plotDur = zeros(1, datSiz(2));
plotPos = zeros(1, datSiz(1));

for jj=1:datSiz(2) %over all durations
    for ii=1:datSiz(1)
        for ww=1:datSiz(3)
            for vv=1:datSiz(4)
                if ~singleBarSt1(ii,jj,ww,vv).empty
                    plotDur(jj) = singleBarSt1(ii,jj,ww,vv).data.table.stimDur;
                    plotPos(ii) = singleBarSt1(ii,jj,ww,vv).data.table.position;
                end
            end
        end
    end
end

plotDurInd = find(ismember(allStimDur, plotDur)); 


% axh = gobjects(datSiz(3), datSiz(1)-1, datSiz(2)-1); % -1 when you add
% the mormalized max and min 
axh = gobjects(datSiz(4), datSiz(3), datSiz(1), datSiz(2));

pCol = cbrewer('qual', 'Set1', 2*(datSiz(4)));


for vv=1:datSiz(4)

    for ww=1:datSiz(3) 

        posCell = generatePositionCell(0.05, 0.975, 0.025, 0.975, 0.005, 0.005, [datSiz(1), datSiz(2)]); % -1 

        axh(vv,ww, :, :) = zeros(datSiz(1), datSiz(2)); %-1

        figure

        for ii=1:datSiz(1) % datSiz-1, since summary stat are in the empty end of the structure

            for jj=1:datSiz(2)

                axh(vv,ww,ii,jj) = axes('position', posCell{ii,jj});
                hold on

                if singleBarSt1(ii,jj,ww,vv).empty
                    continue
                end
                
                plotTime1 = singleBarSt1(ii,jj,ww,vv).subData.baseSub(:,1);
                plotTime2 = singleBarSt2(ii,jj,ww,vv).subData.baseSub(:,1);
                plotDat1 = singleBarSt1(ii,jj,ww,vv).subData.baseSub(:,2);
                plotDat2 = singleBarSt2(ii,jj,ww,vv).subData.baseSub(:,2);
                
                datLen1 = length(plotTime1);
                datLen2 = length(plotTime2);

                plot(plotTime1, plotDat1 , 'linewidth', 2, 'color', pCol(1,:))
                plot(plotTime2, plotDat2, 'linewidth', 2, 'color', pCol(2,:))
                
                if datLen1 <= datLen2
                    diffDat = plotDat1 - plotDat2(1:datLen1);
                    diffTime = plotTime1;
                else
                    diffDat = plotDat1(1:datLen2) - plotDat2;
                    diffTime = plotTime2;
                end
                
                plot(diffTime, diffDat, 'linewidth', 2, 'color', 'k')
                


                line(xRange(jj, :), [0, 0], 'color', [1,1,1]*0.8, 'linewidth', 2)


                if jj==1
                    title(num2str(plotPos(ii)))
                end

                hold off


            end

        end

    end

end


for vv=1:datSiz(4)
    for ww=1:datSiz(3)
        for ii=1:datSiz(2)
            yyLim = vertcat(axh(vv,ww, :,ii).YLim);
            for jj=1:datSiz(1)
                axh(vv,ww,jj,ii).YLim = [min(yyLim(:,1)), max(yyLim(:,2))];
                axh(vv,ww,jj,ii).XLim = xRange(plotDurInd(ii),:);
            end
        end

        set(axh(vv,ww,2:end, :), 'yticklabel', {})
        set(axh(vv,ww,:, 1:end-1), 'xticklabel', {})
    end
end




if nargout==1
    varargout{1} = axh;
end





end