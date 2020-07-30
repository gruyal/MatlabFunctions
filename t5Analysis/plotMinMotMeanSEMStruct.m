function varargout = plotMinMotMeanSEMStruct(minMotMeanSEMStruct)

% function varargout = plotMinMotMeanSEMStruct(minMotMeanSEMStruct)
%
% this function plots the results from allCellMeanMinMotScript (not a
% function). Script is in T5analysis.../minMot...
%
% Function looks for shared FB and SB postions and plots those only


close all

meanDat = minMotMeanSEMStruct.data.mean; 
semDat = minMotMeanSEMStruct.data.sem; 
meanLin = minMotMeanSEMStruct.lin.mean; 
semLin = minMotMeanSEMStruct.lin.sem; 
timeV = minMotMeanSEMStruct.timeVec; 


fbPos = minMotMeanSEMStruct.relFBPos; 
sbPos = minMotMeanSEMStruct.relSBPos; 

sharedPos = intersect(fbPos, sbPos); 

sFBPosInd = zeros(length(sharedPos), 1);
sSBPosInd = zeros(length(sharedPos), 1);

for ii=1:length(sharedPos)
    sFBPosInd(ii) = find(fbPos == sharedPos(ii));
    sSBPosInd(ii) = find(sbPos == sharedPos(ii));
end

dSiz = size(meanDat); 

% posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.02, 0.02, [dSiz(2), dSiz(3)]);
posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.02, 0.02, [length(sharedPos), length(sharedPos)]);
axh = gobjects(size(posCell));

pCol = cbrewer('qual', 'Paired', 6);
patchXX = [1:dSiz(1), dSiz(1):-1:1]';
tPatchXX = timeV(patchXX); 

figure

for ii=1:length(sharedPos) % dSiz(2)
    
    for jj=1:length(sharedPos) %dSiz(3)
        
        axh(ii,jj) = axes('position', posCell{ii,jj});
        hold on 
        
        relFBI = sFBPosInd(ii);
        relSBI = sSBPosInd(jj);
        
        if jj == 1
            title(['FBPos:' num2str(fbPos(relFBI))])
        end
        
        if ii == 1
            axh(ii,jj).YLabel.String = ['SBPos:', num2str(sbPos(relSBI))]; 
        end
        
        tempMD = meanDat(:,relFBI,relSBI); 
        tempSD = semDat(:,relFBI,relSBI); 
        
        tempML = meanLin(:,relFBI,relSBI); 
        tempSL = semLin(:,relFBI,relSBI); 
                
        patchDY = [tempMD + tempSD; flipud(tempMD - tempSD)];
        patchLY = [tempML + tempSL; flipud(tempML - tempSL)];
        
        if sum(tempML(1:100)) ~= 0 % not empty
        
            patch(tPatchXX, patchLY, pCol(1, :), ...
                  'edgecolor', pCol(1, :), 'facealpha', 0.4)
            plot(timeV, tempML, 'linewidth', 2, 'color', pCol(2, :))
            
        end
          
        if sum(tempMD(1:100)) ~=0
        
            axh(ii,jj).Tag = 'used';
            
            patch(tPatchXX, patchDY, pCol(5, :), ...
                  'edgecolor', pCol(5, :), 'facealpha', 0.4)

            plot(timeV, tempMD, 'linewidth', 2, 'color', pCol(6, :))
            
        end
        
        hold off
        
        
    end
    
end


yyMax = [];
yyMin = [];
xxLim = [-250 , 1000]; % a bit wierd since it is downsampled when generated

    
allLinesH = findobj(axh, 'type', 'line', '-and', 'lineWidth', 2);
    
for ii=1:length(allLinesH)
    yDat = allLinesH(ii).YData; 
    yyMax = [yyMax , max(yDat)];
    yyMin = [yyMin , min(yDat)];
end

totMax = max(yyMax); 
totMin = min(yyMin); 

yBuf = (totMax - totMin)/10; 
totYLim = [totMin - yBuf, totMax + yBuf];

for ii=1:length(sharedPos)
    
    for jj=1:length(sharedPos)
        
        axh(ii,jj).XColor = 'none';
        axh(ii,jj).YColor = 'none';

        axh(ii,jj).YLim = totYLim; 
        axh(ii,jj).XLim = xxLim;
        
        if ii==1
            axh(ii,jj).YColor = 'k';
        end
        
        if jj==dSiz(3)
            axh(ii,jj).XColor = 'k';
        end
        
    end
    
end
        






if nargout == 1
    varargout{1} = axh;
end





end

