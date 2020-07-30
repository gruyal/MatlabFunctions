


function varargout = plotMinMotMeanPDND(minMotMeanSEMStruct)

% function varargout = plotMinMotMeanPDND(minMotMeanSEMStruct)
%
% this function plots the results from allCellMeanMinMotScript (not a
% function). Script is in T5analysis.../minMot...
%
% Unlike plotMinMotMeanSEMStruct, this function plots PD and ND from the
% same positions overlayed. The 2 axes in this case are normFBPos and
% nornPosDiff

close all

meanDat = minMotMeanSEMStruct.data.mean; 
semDat = minMotMeanSEMStruct.data.sem; 


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

posCell = generatePositionCell(0.05, 0.975, 0.025, 0.95, 0.02, 0.02, [dSiz(2)-1, dSiz(3)-1]); % -1 since the diagonal is not used
axh = gobjects(size(posCell));

pCol = cbrewer('qual', 'Paired', 6);
patchXX = [1:dSiz(1), dSiz(1):-1:1]';
tPatchXX = timeV(patchXX); 

figure

for ii=1:length(sharedPos)
    
    for jj=ii+1:length(sharedPos)
        
        relFBI = sFBPosInd(ii);
        relSBI = sSBPosInd(jj);
        
        FBV = fbPos(relFBI);
        SBV = sbPos(relSBI);
        
        % to find the reciprocal pair when positions are not identical
        relFBI2 = find(fbPos == SBV);
        relSBI2 = find(sbPos == FBV);

        
        tempMPD = meanDat(:,relFBI,relSBI); 
        tempSPD = semDat(:,relFBI,relSBI); 
        
        tempMND = meanDat(:,relFBI2,relSBI2); 
        tempSND = semDat(:,relFBI2,relSBI2); 
        
                
        patchPDY = [tempMPD + tempSPD; flipud(tempMPD - tempSPD)];
        patchNDY = [tempMND + tempSND; flipud(tempMND - tempSND)];
        
        if sum(tempMND(1:100)) ~= 0 % not empty
            
            axh(ii,jj-1) = axes('position', posCell{ii,jj-1});
            hold on 
            sTit = sprintf('Pos:%d and %d', FBV, SBV);
            title(sTit)
        
            patch(tPatchXX, patchNDY, pCol(1, :), ...
                  'edgecolor', pCol(1, :), 'facealpha', 0.4)
            plot(timeV, tempMND, 'linewidth', 2, 'color', pCol(2, :))
            
            axh(ii, jj-1).XColor = 'none';
            axh(ii, jj-1).YColor = 'none';
        
            axh(ii,jj-1).Tag = 'used';
            
            patch(tPatchXX, patchPDY, pCol(5, :), ...
                  'edgecolor', pCol(5, :), 'facealpha', 0.4)

            plot(timeV, tempMPD, 'linewidth', 2, 'color', pCol(6, :))
            
            hold off
            
        end
                
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

relAx = findobj(axh, 'tag', 'used');

for ii=1:length(relAx)

    
    relAx(ii).YLim = totYLim; 
    relAx(ii).XLim = xxLim;
    
%     for jj=1:dSiz(3)
%         
%         axh(ii,jj).XColor = 'none';
%         axh(ii,jj).YColor = 'none';
% 
%         
%         if ii==1
%             axh(ii,jj).YColor = 'k';
%         end
%         
%         if jj==dSiz(3)
%             axh(ii,jj).XColor = 'k';
%         end
%         
%     end
    
end
        

tempA = findobj(axh(1,:), 'Tag', 'used');
for ii=1:length(tempA)
    tempA(ii).YColor = 'k';
end

tempA = findobj(axh(:, end), 'Tag', 'used');
for ii=1:length(tempA)
    tempA(ii).XColor = 'k';
end





if nargout == 1
    varargout{1} = axh;
end





end

