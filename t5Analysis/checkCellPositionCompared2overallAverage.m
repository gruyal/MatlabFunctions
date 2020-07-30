function varargout = checkCellPositionCompared2overallAverage(cellSBStruct, cellMMStruct, allSBMat, allSBMatTime)

% function varargout = checkCellPositionCompared2overallAverage(cellSBStruct, allSBMat)
%
% This function compares the singlebar data from a single cell to the
% aligned average by plotting the same positions from duration 160ms and width 2 dark bars. 
%
% INPUT
%
% cellSBStruct -            singlebar protocol structure to be fed into generateAlignedSingleBarStwMinMaxDiffWandV
% cellMMStruct -            minimal motion protocol structure. Function
%                           will use only the diagonal positions that are both dark and width 2
% allSBMat -                average of all the singlebar protocols.
%                           Generated in respByPosMeanAndTimeVecScript and saved as T5respByPosMeanAndTimeVec
% allSBMatTime -            time vactor for allSBMat. Generate in the same
%                           script
%
% OUTPUT
%
% axh -                     optional. Axes handles to the plots

close all

singleBarSt = generateAlignedSingleBarStwMinMaxDiffWandV(cellSBStruct);

rDur = [0.04, 0.16]';
rWid = 2;
rVal = 0;
%single bar data
relTab = cellSBStruct.gratingTable; 
uPos = unique(relTab.position);
uDur = unique(relTab.stimDur);
uVal = unique(relTab.value);
uWid = unique(relTab.width);

relSBSt = singleBarSt(1:end-1, ismember(uDur, rDur), uWid == rWid, uVal == rVal); 

relParams = singleBarSt(end, end, end, end);
relMaxExt = relParams.maxExtPosVal;
relPD = sign(relParams.maxInhPosVal - relMaxExt);

preRelPos = (uPos - relMaxExt); 
if relPD == 1
    relPos = preRelPos * relPD;
else
    relPos = (preRelPos - rWid +1) * relPD; 
end

% minMot data

mmRes = calcMinMotExtLinCompDiffWandVwTable(cellMMStruct, cellSBStruct);
mmTab = mmRes.mmTable; 

tabInd1 = mmTab.FBPos == mmTab.SBPos;
tabInd2 = mmTab.width == rWid;
tabInd3 = mmTab.FBVal == rVal & mmTab.SBVal == rVal;
tabIndS = mmTab.timeDiff == rDur(1);
tabIndL = mmTab.timeDiff == rDur(2);

totSI = tabInd1 & tabInd2 & tabInd3 & tabIndS; 
totLI = tabInd1 & tabInd2 & tabInd3 & tabIndL; 

mmShortData = mmRes.mmResult(totSI); 
mmLongData = mmRes.mmResult(totLI); 

mmSPos = mmTab.normFBPos(totSI); 
mmLPos = mmTab.normFBPos(totLI); 

% all single bar average
sbTDOrd = [0.02, 0.04, 0.08, 0.16, 0.32];
sbWidOrd = [1,2,4];
sbValOrd = [0,1];
relAllMeanLongD = allSBMat(:,:,sbTDOrd == rDur(2), sbWidOrd == rWid, sbValOrd == rVal); % 160ms, wid 2 , val 0
relAllMeanShortD = allSBMat(:,:,sbTDOrd == rDur(1), sbWidOrd == rWid, sbValOrd == rVal); % 40ms, wid 2 , val 0
relLen = size(relAllMeanLongD,1);
cenPos = floor(size(relAllMeanLongD,2)/2);
relAllPosInd = find(sum(isnan(relAllMeanLongD)) < relLen);
cenRelAllInd = relAllPosInd - cenPos; 

posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.005, 0.02, [length(relAllPosInd), 2]);
axh = gobjects(size(posCell));

figure('position', [200, 650, 2300, 600])

pCol = cbrewer('qual', 'Set1', 4);



for ii=1:length(relAllPosInd)
    
    
    % plotting 40ms 
    axh(ii, 1) = axes('position', posCell{ii,1});
    hold on 
    tempAllDatS = relAllMeanShortD(:, relAllPosInd(ii)); 
    
    plot(allSBMatTime, tempAllDatS, 'linewidth', 2, 'color', 'k')
    
    relI = relPos == cenRelAllInd(ii); 
    if sum(relI) ~= 0 && ~isempty(relSBSt) && ~relSBSt(relI,1).empty
        relCellDat = relSBSt(relI,1).subData.baseSub;
        plot(relCellDat(:,1), relCellDat(:,2), 'linewidth', 2, 'color', pCol(1,:))
    end
    
    relMMIS = mmSPos == cenRelAllInd(ii); 
    if sum(relMMIS) == 1
        mmDatS = mmShortData(relMMIS).subData.baseSub;
        plot(mmDatS(:,1), mmDatS(:,2), 'linewidth', 2, 'color', pCol(2,:))
    end
    
    title({'normPos:'; num2str(cenRelAllInd(ii))})
    hold off
    
    % plotting 160ms
    axh(ii, 2) = axes('position', posCell{ii,2});
    hold on 
    tempAllDatL = relAllMeanLongD(:, relAllPosInd(ii)); 
    plot(allSBMatTime, tempAllDatL, 'linewidth', 2, 'color', 'k')

    if sum(relI) ~= 0 && ~isempty(relSBSt) && ~relSBSt(relI,2).empty
        relCellDat = relSBSt(relI,2).subData.baseSub;
        plot(relCellDat(:,1), relCellDat(:,2), 'linewidth', 2, 'color', pCol(1,:))
    end
    
    relMMIL = mmLPos == cenRelAllInd(ii); 
    if sum(relMMIL) == 1
        mmDatL = mmLongData(relMMIL).subData.baseSub;
        plot(mmDatL(:,1), mmDatL(:,2), 'linewidth', 2, 'color', pCol(2,:))
    end
    
    hold off
    
    
end


% adding legend
yyLeg = 5;
xxLeg = 100;

set(gcf, 'currentAxes', axh(1,2))

text(xxLeg, yyLeg, 'allAverage','color', 'k', 'fontSize', 12)
text(xxLeg, yyLeg + 4, 'singleBar','color', pCol(1,:), 'fontSize', 12)
text(xxLeg, yyLeg + 8, 'minMot','color', pCol(2,:), 'fontSize', 12)


%orginizing scales
xxRange = [-200, 500; -200, 750];

yyMax = cell(1,2);
yyMin = cell(1,2);

for ii=1:2
    
    relLineH = findobj(axh(:,ii), 'Type', 'Line'); 
    for jj=1:length(relLineH)
        
        yDat = relLineH(jj).YData; 
        yyMax{ii} = [yyMax{ii} , max(yDat)];
        yyMin{ii} = [yyMin{ii} , min(yDat)];
        
    end
    
end
    
totMax = cellfun(@max, yyMax); 
totMin = cellfun(@min, yyMin); 

rMax = ceil(totMax/10)*10; 
rMin = floor(totMin/10)*10; 

rRange = rMax - rMin;
yyRange = [rMin - rRange/20; rMax + rRange/10]'; 

for ii=1:size(axh,2)
    
    for jj=1:size(axh,1)
        
        
        axh(jj,ii).XLim = xxRange(ii, :);
        axh(jj,ii).YLim = yyRange(ii,:);
        
        if jj > 1
            axh(jj,ii).YColor = 'none'; 
        end
        
%         if ii==1
%             axh(jj,ii).XColor = 'none'; 
%         end
        
    end
    
end


if nargout == 1
    varargout{1} = axh;
end




end