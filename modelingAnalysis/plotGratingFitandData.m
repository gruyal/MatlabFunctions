function varargout  = plotGratingFitandData(allCellDat)

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
% 


relDur = [40, 160]; % since this is common for all T5 cells 

[mgSt, mgFitTab] = fitCosGratingDatAndModel(allCellDat);

snL = unique(mgFitTab.index(mgFitTab.direction == 1 & ismember(mgFitTab.duration, relDur))); 

    

figure('position', [1600, 200, 825, 1000])


posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.02, 0.05, [3, length(snL)]);
axh = gobjects(size(posCell)); 
yyLims = zeros([size(axh), 2]); 

for si=1:length(snL)

    sn = snL(si); 

    tempTab = mgFitTab(mgFitTab.index == sn & mgFitTab.dataFlag == 1, :);
    tempDir = [0, 1];
    tempDur = tempTab.duration; 
    tempStP = tempTab.phase; 

    for ii=1:3
        axh(ii, si) = axes('position', posCell{ii, si}); 
        hold on 
    end

    preCol = cbrewer('qual', 'Paired', 6); 
    pCol = preCol([1,2,5,6], :);
    dColInd = [2,4];
    mColInd = [1,3];

    for dd=1:length(tempDir)

        tempT = mgFitTab(  mgFitTab.dataFlag == 1 & ...
                           mgFitTab.direction == tempDir(dd) & ...
                           mgFitTab.duration == tempDur & ...
                           mgFitTab.phase == tempStP, :);
                       
        tempInd = tempT.index; 

        relT = allCellDat.MG(tempInd).time;
        relD = allCellDat.MG(tempInd).data; 
        relM = allCellDat.MG(tempInd).model;

        datFit = mgSt.MG(tempInd).dataFit;
        modFit = mgSt.MG(tempInd).modelFit;

        set(gcf, 'currentaxes', axh(1,si))
        plot(relT, relD, 'color', pCol(dColInd(dd),:), 'linewidth', 2)
        plot(datFit)
        fitLH = axh(1, si).Children(1);  
        fitLH.Color = pCol(dColInd(dd),:); 
        legend('off')

        set(gcf, 'currentaxes', axh(2, si))
        plot(relT, relM, 'color', pCol(mColInd(dd),:), 'linewidth', 2)
        plot(modFit)
        fitLH = axh(2, si).Children(1);  
        fitLH.Color = pCol(mColInd(dd),:); 
        legend('off')

        set(gcf, 'currentaxes', axh(3, si))
        plot(relT, relD, 'color', pCol(dColInd(dd),:), 'linewidth', 2)
        plot(relT, relM, 'color', pCol(mColInd(dd),:), 'linewidth', 2)

    end

    if tempDur == 40
        xLMax = 1200; 
    else
        xLMax = 4000; 
    end

    for ii=1:3
        axh(ii, si).XLim = [-100, xLMax];
        axh(ii, si).XLabel.Visible = 'off';
        axh(ii, si).YLabel.Visible = 'off';
        hold(axh(ii, si), 'off')
        yyLims(ii, si, :) = axh(ii, si).YLim; 
    end
    clear temp* rel* datFit modFit timeFit

end

axS = size(axh);
titText = {'Data v Fit'; 'Model v Fit'; 'Data v Model'};
titText2 = {'stepDur 40ms'; 'stepDur 160ms'};

yMins = min(yyLims(:,:,1));
yMaxs = max(yyLims(:,:,2));

for ii=1:axS(1)
    for jj=1:axS(2)

        axh(ii,jj).YLim = [yMins(jj), yMaxs(jj)];

        if ii > 1
            axh(ii,jj).YColor = 'none';
        end

        if mod(jj,2)
            axh(ii,jj).XColor = 'none';
        end

        if jj == 1
            axh(ii,jj).Title.String = titText{ii};
        end

        if ismember(jj, [2,4]) && ii==2
            axh(ii,jj).Title.String = titText2{jj/2};
            titPos = axh(ii,jj).Title.Position;
            axh(ii,jj).Title.Position = titPos + [0, titPos(2)/10, 0];
        end

    end
end


lH1 = findobj(axh(end,end), 'color', pCol(2,:));
lH2 = findobj(axh(end,end), 'color', pCol(4,:));

legend([lH1, lH2], 'ND', 'PD')
legend('boxoff')

if nargout == 1
    varargout{1} = axh;
end

    
end

    
