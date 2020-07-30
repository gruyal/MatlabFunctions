function varargout  = plotGratingFitandDataFromAlignedDat(alignGratSt)

% function varargout  = plotGratingFitandData(allCellDat)
%
% This function plots grating data and model after it fit each with the
% optimal cosine
%
% INPUT
% allCellDat -      generated using orginizingClusterData function 
%
% OUTPUT
% axh -             optional, handles to all the axes 
% mgFitTab -        optional, table with fit results


[mgSt, mgFitTab] = fitCosAlignedGratingSt(alignGratSt);
    
relDir = [-1, 1];

preCol = cbrewer('qual', 'Paired', 6); 
pCol = preCol([1,2,5,6], :);
dColInd = [2,4];

xxRange = [-250, 2000; -250, 7000];

uComb = unique(mgFitTab.grtComb); 
uDur = unique(mgFitTab.stimDur); 
uStP = unique(mgFitTab.startPhase); 

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.02, 0.05, [length(uStP), length(uDur)]);
axh = gobjects(length(uComb), length(uDur), length(uStP)); 
yyLims = zeros([size(axh), 2]); 

for comb=1:length(uComb)

    figure('position', [1600, 200, 825, 1000], 'name', ['Comb:', num2str(comb)])
    
    tempTab = mgFitTab(mgFitTab.grtComb == uComb(comb), :);
    
    
    for dd=1:length(uDur)
        
        for st=1:length(uStP)
            
            axh(comb, dd, st)  = axes('position', posCell{st, dd});
            hold on 
            
            for rr=1:length(relDir)

                
                tempT = tempTab(  tempTab.direction == relDir(rr) & ...
                                  tempTab.stimDur == uDur(dd) & ...
                                  tempTab.startPhase == uStP(st), :);

                tempInd = tempT.index; 

                relT = alignGratSt.result(tempInd).subData.baseSub(:,1);
                relD = alignGratSt.result(tempInd).subData.baseSub(:,2);
            

                datFit = mgSt(tempInd).Fit;
                fitTime  = mgSt(tempInd).RelTime;
                
                plot(fitTime, zeros(size(fitTime)), 'color', [1,1,1]*0.8, 'linewidth', 4)
                plot(relT, relD, 'color', pCol(dColInd(rr)-1,:), 'linewidth', 1)
                plot(datFit)
                fitLH = axh(comb, dd, st).Children(1);  
                fitLH.Color = pCol(dColInd(rr),:); 
                fitLH.LineWidth = 2; 
                legend('off')
                
            end
            
            hold off
            
            yyLims(comb, dd, st, :) = axh(comb, dd, st).YLim; 
            
        end
                
    end
    
end
                        

axS = size(axh);        

yMins = min(yyLims(:,:,:,1));
yMaxs = max(yyLims(:,:,:,2));

for ii=1:axS(1)
    
    for jj=1:axS(2)
        
        for kk=1:axS(3) 

            axh(ii,jj,kk).YLim = [yMins(kk), yMaxs(kk)];
            axh(ii,jj,kk).XLim = xxRange(jj, :); 

            if kk > 1
                axh(ii,jj,kk).YColor = 'none';
            end

            if jj == 1
                axh(ii,jj,kk).Title.String = ['startPhase:', num2str(uStP(kk))];
            end

            if kk == 1
                axh(ii,jj,kk).YLabel.String = ['StepDur:', num2str(uDur(jj))];
            end
            
        end

    end
    
    lH1 = findobj(axh(comb, end,end), 'color', pCol(2,:));
    lH2 = findobj(axh(comb, end,end), 'color', pCol(4,:));

    legend([lH1(end), lH2(end)], 'ND', 'PD')
    legend('boxoff')
    
end







if nargout == 1
    varargout{1} = axh;
elseif nargout == 2
    varargout{1} = axh;
    varargout{2} = mgFitTab;
end
    
end

    
