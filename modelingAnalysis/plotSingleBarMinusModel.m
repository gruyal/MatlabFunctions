function varargout = plotSingleBarMinusModel(allCellDat)

% function varargout = plotSingleBarMinusModel(allCellDat)
%
% This function plots all the single bar data as a heatmap by position,
% width and duration , showing the difference between data and model to see
% if there are any consistent biases
%
% INPUT
% allCellDat -      generated using organizingClusterData
%
% OUTPUT
% axh -             (optional) handle to the axes in the figure




relWid = [1,2,4];
relDur = [40, 160];
relPos = -6:6; 
xxTick = 55:50:205;

numCol = 11; 

xSt1 = 0.075; 
xEnd1 = 0.975; 
ySt1 = 0.1; 
yEnd1 = 0.95; 

posCell = generatePositionCell(xSt1, xEnd1, ySt1, yEnd1, 0.02, 0.02, [length(relDur), length(relWid)]);
axh = gobjects(size(posCell)); 

cMap = flipud(cbrewer('div', 'RdYlBu', numCol));

sbTable = allCellDat.table(allCellDat.table.protType == 0, :); 

figure('position', [1750, 250, 600, 1100])
colormap(cMap)
allMax = zeros(size(posCell)); 
allMin = zeros(size(posCell)); 

for ww=1:length(relWid)

    for dd=1:length(relDur)

        tempLen = sbTable.dataLength(sbTable.width == relWid(ww) & sbTable.duration == relDur(dd));  
        maxLen = max(tempLen); 

        preImage = nan(length(relPos), maxLen); 

        for pp=1:length(relPos)

            relInd = sbTable.index(sbTable.width == relWid(ww) & sbTable.duration == relDur(dd) & sbTable.position == relPos(pp));  
            relLen = sbTable.dataLength(sbTable.width == relWid(ww) & sbTable.duration == relDur(dd) & sbTable.position == relPos(pp));  

            if isempty(relInd)
                continue
            end

            assert(length(relInd) == 1, 'index wrong')

            dataV = allCellDat.SB(relInd).data;
            modelV = allCellDat.SB(relInd).model;
            timeV = allCellDat.SB(relInd).time;

            preImage( pp , 1:relLen) = dataV - modelV; 

        end

        if isempty(preImage)
            continue
        end

        axh(dd,ww) = axes('position', posCell{dd,ww});

        imAlpha = ones(size(preImage));
        imAlpha(isnan(preImage)) = 0;

        imagesc(preImage, 'AlphaData',imAlpha)
        line([0, maxLen], [1,1]*length(relPos)/2 + 0.5, 'color', [1,1,1]*0.5, 'linewidth', 2)

        allMax(dd,ww) = nanmax(preImage(:));
        allMin(dd,ww) = nanmin(preImage(:));



    end

end

relxxTick = xxTick(xxTick <= length(timeV)); 

axS = size(axh); 
for ii=1:axS(1)
    for jj=1:axS(2)

        if isempty(findobj(axh(ii,jj), 'type', 'Axes')) % in case there is a graphics place holder there
            continue
        end

        axh(ii,jj).CLim = [min(allMin(:)), max(allMax(:))]; 
        axh(ii,jj).Color = [1,1,1]*0.85;
        axh(ii,jj).YTick = 1:2:length(relPos); 
        axh(ii,jj).YTickLabel = arrayfun(@num2str, relPos(1:2:end), 'uniformoutput', 0);
        axh(ii,jj).XTick = relxxTick; 
        axh(ii,jj).FontSize = 14; 
        axh(ii,jj).XTickLabel = arrayfun(@(x) num2str(floor(x)), timeV(relxxTick), 'uniformoutput', 0); 

        if ii > 1
            axh(ii,jj).YColor = 'none';
        end

        if jj < axS(2)
            axh(ii,jj).XColor = 'none';
        end

        if ii==1
            axh(ii,jj).YLabel.String = ['Wid:', num2str(relWid(jj))];
            axh(ii,jj).YLabel.FontWeight = 'bold';
        end

        if jj==1
            axh(ii,jj).Title.String = ['Dur:', num2str(relDur(ii))];
        end

    end
end

axh2 = axes('position', [xSt1, ySt1 - 0.065, xEnd1-xSt1, 0.03]);
colBarDat = linspace(min(allMin(:)), max(allMax(:)), numCol); 
xTick2 = [1,3,6,9,11];  % adjust to numCol 

imagesc(colBarDat)
axh2.YColor = 'none'; 
axh2.XTick = xTick2;
axh2.FontSize = 14; 
axh2.XTickLabel = arrayfun(@(x) num2str(x, 2), colBarDat(xTick2), 'uniformoutput', 0); 



if nargout == 1
    varargout{1} = axh; 
end




end
            
            
            
            