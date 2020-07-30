function plotSlopeFitResultsOnData(sbStruct, fitTab, parSt)

% function plotSlopeFitResultsOnData(sbStruct, fitTab, parSt)
%
% This function is specifically designed to check results from
% 'checkingDiff4WStatScript'. it takes a specific presentation and plots
% the corresponding fit on it
%
% INPUT
% sbStruct -        singleBarSt created in checkingDiff4WStatScript
% fitTab-           tausTab created in the same script
% parSt -           structure that provides details on which stim response
%                   to plot. W4 and Dur160 are assumed. cellNum, val and
%                   nPos (normalized position - as in tausTab) need to be specified 

cellNum = parSt.cellNum;
val = parSt.val;
nPos = parSt.nPos;

relRow = fitTab(fitTab.cellNum == cellNum & fitTab.nPos == nPos & fitTab.val == val, :);

if isempty(relRow)
    warning('specific combination is empty')
elseif height(relRow) > 1
    warning('more than one row for specific combination')
end

relWid = 4; 
relDur = 0.16;

tempStat = sbStruct(cellNum).result(end,end,end,end);
tempPos = sbStruct(cellNum).positions;
   
tempEPos = tempStat.maxExtPosVal;
tempIPos = tempStat.minInhPosVal;
normPos = nan(size(tempPos)); 

relPD = sign(tempIPos - tempEPos);

if relPD == 1
    normPos = tempPos - tempEPos;
else
    normPos = (tempPos - tempEPos - relWid + 1) * relPD; 
end 

dInd = sbStruct(cellNum).durations == relDur;
wInd = sbStruct(cellNum).widths == relWid; 
vInd = sbStruct(cellNum).vals == val;
pInd = normPos == nPos;  

tempTime = sbStruct(cellNum).result(pInd, dInd, wInd, vInd).subData.baseSub(:,1);
tempDat = sbStruct(cellNum).result(pInd, dInd, wInd, vInd).subData.baseSub(:,2);


riseTime = tempTime(relRow.riseStInd:relRow.riseEndInd);
decayTime = tempTime(relRow.decayStInd:relRow.decayEndInd);
rSlope = relRow.riseSlope;
rInt = relRow.riseInt; 
dSlope = relRow.decaySlope;
dInt = relRow.decayInt; 

clf

hold on 

plot(tempTime, tempDat, 'linewidth', 1, 'color', [1,1,1]*0.7)
plot(riseTime, rSlope * riseTime + rInt, 'linewidth', 2, 'color', 'r')
plot(decayTime, dSlope * decayTime + dInt, 'linewidth', 2, 'color', 'b')

hold off

title(['c:', num2str(cellNum), ' v:', num2str(val), ' np:', num2str(nPos)])


end