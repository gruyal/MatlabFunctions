function varargout = plotAlignedFlightDataOL(alignExpStruct, stimNum)

% function varargout = plotAlignedFlightDataOL(alignExpStruct)
% 
% This function uses the output from alignFlightArenaExpStruct to plot the
% data from stimulus stimNum. Plot include recalculated L-R, L+R and Freq
% data for individual trials and mean. 
%
% 


close all

timeCh = 1; lCh = 3; rCh = 4; freqCh = 5;  % ch 2 is L-R from the flight controller (is used for closed loop but less accurate)

relData = alignExpStruct.dataOL(stimNum, end);

excludedT = ~relData.usefulFlag;

datTime = relData.alignData(:,:,timeCh);
lMinusR = relData.alignData(:,:, lCh) - relData.alignData(:,:, rCh);
lPlusR = relData.alignData(:,:, lCh) + relData.alignData(:,:, rCh);
freqDat = relData.alignData(:,:, freqCh); 


meanTime = relData.meanData(:,timeCh);
meanLminR = relData.meanData(:,lCh) - relData.meanData(:, rCh);
meanLplusR = relData.meanData(:,lCh) + relData.meanData(:, rCh);
meanFreq = relData.meanData(:,freqCh);

exFlag = 0;
if sum(excludedT > 0)
    exTime = datTime(:, excludedT);
    exLMR = lMinusR(:,excludedT);
    exLPR = lPlusR(:,excludedT);
    exFreq = freqDat(:, excludedT);
    exFlag = 1;
end

posData = relData.expPos;

rCols = cbrewer('qual', 'Paired', 8);

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -0.1, 0.02, 3); 
axh = gobjects(1,3);

% L-R

axh(3) = axes('position', posCell{1});
hold on 
plot(datTime, lMinusR, 'linewidth', 2, 'color', rCols(1, :))
title('L-R')

yyLim = axh(3).YLim;
yyLimNew = yyLim(1) - 2; % shifts it below data

plot(posData(:,1), posData(:,2)./10 + yyLimNew, 'color', [1,1,1]*0.7)
plot(meanTime, meanLminR, 'linewidth', 2, 'color', rCols(2, :))

if exFlag
    plot(exTime, exLMR, 'linewidth', 1, 'color', [1,1,1]*0.6);
end

hold off

axh(3).YLim = [yyLimNew, yyLim(2)];
% L+R
axh(2) = axes('position', posCell{2});
hold on 
plot(datTime, lPlusR, 'linewidth', 2, 'color', rCols(3, :))
plot(meanTime, meanLplusR, 'linewidth', 2, 'color', rCols(4, :))

if exFlag
    plot(exTime, exLPR, 'linewidth', 1, 'color', [1,1,1]*0.6);
end

title('L+R')
hold off


axh(1) = axes('position', posCell{3});
hold on 
plot(datTime, freqDat, 'linewidth', 2, 'color', rCols(5, :))
plot(meanTime, meanFreq, 'linewidth', 2, 'color', rCols(6, :))

if exFlag
    plot(exTime, exFreq, 'linewidth', 1, 'color', [1,1,1]*0.6);
end


title('Frequency')
hold off


if nargout > 0
    varargout{1} = axh;
end


end



