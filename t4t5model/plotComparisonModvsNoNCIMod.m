function plotComparisonModvsNoNCIMod(resSt, resTab, plotPar)

% function plotComparisonModvsNoNCIMod(resSt, resTab, plotPar)
%
% this function compares results from removingNCIScriptV2. it compares the
% simpleV2 model with the same parameters after non-preferred I is removed
% (zeroed). 
%
% INPUT
%
% resSt -           resultStucture from removingNCIScriptV2 (contains V, ge, gi for each
%                   simulation)
% resTab -          resultsTable from same file. contains:  count    cellT    cellN    relIter    barVal    barWid    stepD    PDFlag    fullModFlag    optResVal
%                   count in table refers to the position in the structure              
% plotPar -         plotting parameters (same varaibles as in table) to
%                   select which data to plot from resSt. 
%                   sould contain the following fields
%                   plotPar.cellN, plotPar.cellT, plotPar.iterOrd,
%                   plotPar.barVal, plotPar.barWid
%


allSteps = unique(resTab.stepD); 
numSteps = length(allSteps);

movMeanWin = 51; 

relInds = resTab.cellN == plotPar.cellN & resTab.cellT == plotPar.cellT & resTab.iterOrd == plotPar.iterOrd & ...
          resTab.barVal == plotPar.barVal & resTab.barWid == plotPar.barWid; 

tempTab = resTab(relInds, :); 

assert(height(tempTab) == (numSteps * 2 * 2), 'plotPar did not limit data enough')

prePosCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.05, 0.05, [ceil(numSteps/2), 2]);
axh = gobjects([numel(prePosCell), 3, 2]);

pnCol = cbrewer('qual', 'Set1', 3);
eiCol = cbrewer('div', 'PiYG', 4);

figName = sprintf('CN:%d iterO:%d Val:%d Wid:%d', plotPar.cellN, plotPar.iterOrd, plotPar.barVal, plotPar.barWid);
figure('position', [150, 100, 1500, 900], 'name', figName);

for dd=1:numSteps
    
    xPosSt = prePosCell{dd}(1);
    yPosSt = prePosCell{dd}(2);
    xPosEnd = sum(prePosCell{dd}([1,3]));
    yPosEnd = sum(prePosCell{dd}([2,4]));
    
    stepTab = tempTab(tempTab.stepD == allSteps(dd), :);
    
    pdFullInd = stepTab.count(stepTab.PDFlag == 1 & stepTab.fullModFlag == 1);
    ndFullInd = stepTab.count(stepTab.PDFlag == 0 & stepTab.fullModFlag == 1);
    pdNCIInd = stepTab.count(stepTab.PDFlag == 1 & stepTab.fullModFlag == 0);
    ndNCIInd = stepTab.count(stepTab.PDFlag == 0 & stepTab.fullModFlag == 0);
    
    timeV = resSt(pdFullInd).tVec; 
    pdVF = resSt(pdFullInd).vVec;
    pdGEF = resSt(pdFullInd).geVec;
    pdGIF = resSt(pdFullInd).giVec;
    ndVF = resSt(ndFullInd).vVec;
    ndGEF = resSt(ndFullInd).geVec;
    ndGIF = resSt(ndFullInd).giVec;
    pdVN = resSt(pdNCIInd).vVec;
    pdGEN = resSt(pdNCIInd).geVec;
    pdGIN = resSt(pdNCIInd).giVec;
    ndVN = resSt(ndNCIInd).vVec;
    ndGEN = resSt(ndNCIInd).geVec;
    ndGIN = resSt(ndNCIInd).giVec;
    
    [~, maxPDFI] = max(movmean(pdVF, movMeanWin));
    [~, maxNDFI] = max(movmean(ndVF, movMeanWin));
    [~, maxPDNI] = max(movmean(pdVN, movMeanWin));
    [~, maxNDNI] = max(movmean(ndVN, movMeanWin));
    
    tempPosCell = generatePositionCell(xPosSt, xPosEnd, yPosSt, yPosEnd, 0.01, 0.01, [3, 2]);
    
    axh(dd, 1, 1) = axes('position', tempPosCell{1,1});
    hold on 
    plot(timeV, pdVF, 'color', pnCol(1,:), 'linewidth', 2)
    plot(timeV, ndVF, 'color', pnCol(2,:), 'linewidth', 2)
    plot(timeV(maxPDFI), pdVF(maxPDFI), 'color', pnCol(1,:), 'marker', 'o', ...
         'markerfacecolor',  pnCol(1,:), 'markersize', 10)
    plot(timeV(maxNDFI), ndVF(maxNDFI), 'color', pnCol(2,:), 'marker', 'o', ...
         'markerfacecolor',  pnCol(2,:), 'markersize', 10)
    hold off
    
    axh(dd, 1, 2) = axes('position', tempPosCell{1,2});
    hold on 
    plot(timeV, pdVN, 'color', pnCol(1,:), 'linewidth', 2)
    plot(timeV, ndVN, 'color', pnCol(2,:), 'linewidth', 2)
    plot(timeV(maxPDNI), pdVN(maxPDNI), 'color', pnCol(1,:), 'marker', 'o', ...
         'markerfacecolor',  pnCol(1,:), 'markersize', 10)
    plot(timeV(maxNDNI), ndVN(maxNDNI), 'color', pnCol(2,:), 'marker', 'o', ...
         'markerfacecolor',  pnCol(2,:), 'markersize', 10)
    hold off
    
    axh(dd, 2, 1) = axes('position', tempPosCell{2,1});
    hold on 
    plot(timeV, pdGEF, 'color', eiCol(4,:), 'linewidth', 2)
    plot(timeV, pdGIF, 'color', eiCol(1,:), 'linewidth', 2)
    hold off
    title(['StepDur:', num2str(allSteps(dd))])
    
    axh(dd, 3, 1) = axes('position', tempPosCell{3,1});
    hold on 
    plot(timeV, ndGEF, 'color', eiCol(4,:), 'linewidth', 2)
    plot(timeV, ndGIF, 'color', eiCol(1,:), 'linewidth', 2)
    hold off
    
    axh(dd, 2, 2) = axes('position', tempPosCell{2,2});
    hold on 
    plot(timeV, pdGEN, 'color', eiCol(4,:), 'linewidth', 2)
    plot(timeV, pdGIN, 'color', eiCol(1,:), 'linewidth', 2)
    hold off
    
    axh(dd, 3, 2) = axes('position', tempPosCell{3,2});
    hold on 
    plot(timeV, ndGEN, 'color', eiCol(4,:), 'linewidth', 2)
    plot(timeV, ndGIN, 'color', eiCol(1,:), 'linewidth', 2)
    hold off
    
    


end


axSz = size(axh);

for ii=1:axSz(1)
    
    for jj=1:axSz(2)
        
        tempYMax = 0;
        tempYMin = 0;
        for kk=1:axSz(3)
            
            axh(ii,jj,kk).XLim = [0, 1250 * ii + floor(ii/6) * 2000];
            
            if ii==1 && kk==1
                if jj==1
                    legend(axh(ii,jj,kk), 'vPD', 'vND')
                elseif jj==2
                    legend(axh(ii,jj,kk), 'GE', 'GI')
                    
                end
            end
            
            if axh(ii,jj,kk).YLim(2) > tempYMax
                tempYMax = axh(ii,jj,kk).YLim(2);
            end
            
            if jj == 1 && axh(ii,jj,kk).YLim(1) < tempYMin
                tempYMin = axh(ii,jj,kk).YLim(1);
            elseif jj > 1
                tempYMin = -0.1; 
            end
            
            if kk == 1
                axh(ii,jj,kk).XColor = 'none';
            else
                if jj==2
                    axh(ii,jj,kk).XColor = pnCol(1,:);
                elseif jj==3
                    axh(ii,jj,kk).XColor = pnCol(2,:);
                end
                
            end
            
            if ismember(ii, [1,4]) && jj==1
                if kk==1
                    axh(ii,jj,kk).YLabel.String = 'Full';
                elseif kk==2
                    axh(ii,jj,kk).YLabel.String = 'no NC I';
                end
            end
            
        end
        
        axh(ii,jj,1).YLim  = [tempYMin, tempYMax];
        axh(ii,jj,2).YLim  = [tempYMin, tempYMax];
        
    end
    
end







end