function varargout = plotCombSingleBarStFromMinMotAnalysis(minMotAnaSt)

% function axh = plotCombSingleBarStFromMinMotAnalysis(minMotAnaSt)
%
% This function is a QC function to see how singleBar data has been
% combined with the diagonal responses in the minMot protocol. It is
% designed to work with the output from calcMinMotExtLinCompDiffWandVwTable
% that generate a combSingBarSt structure. 
%
% INPUT
% 
% combSingBarSt -       output from calcMinMotExtLinCompDiffWandVwTable
%
% OUTPUT
%
% axh -                 axes handles for all the axes




assert(isfield(minMotAnaSt, 'combSingBarSt'), 'function designed specifically for calcMinMotExtLinCompDiffWandVwTable output')


% plotting combined singleBar and MinMot diagonal

combDat = minMotAnaSt.combSingBarSt; 
allComb = length(combDat); 


singDat = minMotAnaSt.justSingBarSt; 
if isfield(singDat, 'singleBarVec')
    allSing = length(singDat); 
else
    allSing = 0;
end

maxAll = max(allSing, allComb);  % so as to maintain the same size for the axes in the 2 figures

axh = gobjects(1, allComb); 
axNum = ceil(sqrt(maxAll));
posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.02, 0.02, [axNum, axNum]);
pCol = cbrewer('qual', 'Set1', 5);


figure

for ii=1:allComb
    
    axh(ii) = axes('position', posCell{ii});
    
    hold on 
    plot(combDat(ii).singleBarVec, 'linewidth', 2, 'color', pCol(1,:))
    plot(combDat(ii).mmDiagVec, 'linewidth', 2, 'color', pCol(2,:))
    plot(combDat(ii).combined, 'linewidth', 2, 'color', pCol(4,:))
    
    sbTabI = combDat(ii).sbTabIndex;
    tempSBTab = minMotAnaSt.sbTable(minMotAnaSt.sbTable.index == sbTabI, :);
    titString = sprintf('Dur:%.2f Pos:%d Wid:%d Val:%d', ...
                        tempSBTab.stimDur, tempSBTab.position, tempSBTab.width, tempSBTab.value);
    
    axh(ii).XLim = [0.5, 3] * 10^4; 
    title(titString)
    hold off
    
end 


for ii=1:allComb
    axh(ii).XColor = 'none';
    axh(ii).YColor = 'none'; 
end

legend('SB', 'MM', 'MEAN', 'location', 'southwest')

for ii=1:axNum:allComb
    axh(ii).YColor = 'k'; 
end
   
% plotting just the single bar data

axh2 = gobjects(1, allSing); 

if isfield(singDat, 'singleBarVec') % meaning if it isn't empty (since length would give one)
    
    figure

    for ii=1:allSing

        axh2(ii) = axes('position', posCell{ii});

        hold on 
        plot(singDat(ii).singleBarVec, 'linewidth', 2, 'color', pCol(1,:))

        sbTabI = singDat(ii).sbTabIndex;
        tempSBTab = minMotAnaSt.sbTable(minMotAnaSt.sbTable.index == sbTabI, :);
        titString = sprintf('Dur:%.2f Pos:%d Wid:%d Val:%d', ...
                            tempSBTab.stimDur, tempSBTab.position, tempSBTab.width, tempSBTab.value);

        axh2(ii).XLim = [0.5, 3] * 10^4; 
        title(titString)
        hold off

    end  
    
end

for ii=1:allSing
    axh2(ii).XColor = 'none';
    axh2(ii).YColor = 'none'; 
end

for ii=1:axNum:allSing
    axh2(ii).YColor = 'k'; 
end


allAx = [axh, axh2];

allLinesH = findobj(allAx, 'Type', 'Line', '-and', 'lineWidth', 2);
yyMax = zeros(1,length(allLinesH)); 
yyMin = zeros(1,length(allLinesH)); 

for ii=1:length(allLinesH)
    yDat = allLinesH(ii).YData; 
    yyMax(ii) = max(yDat);
    yyMin(ii) = min(yDat);
end

totMax = max(yyMax); 
totMin = min(yyMin); 

yBuf = (totMax - totMin)/10; 
totYLim = [totMin - yBuf, totMax + yBuf];

for ax=1:length(allAx)
    
    allAx(ax).YLim = totYLim;
    
end



if nargout == 1
    varargout{1} = allAx;
end







end

