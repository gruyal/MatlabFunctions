function [gratingFitSt, gratingFitTab] = fitHilbertAlignedGratingSt(alignMovGratSt)

%%% !!! Function not finished - finds false positive in noisy conditions 

% function [gratingFitSt, gratingFitTab] = fitHilbertAlignedGratingSt(pStruct)
%
% This function was created since for some of the cells cosine didn't quite
% capture the shape of the response and therefore the amplitude was
% underestimated. It perform a hilbert transform and uses the zeros and Pi
% positions in the phase to estimate the max and min in the data.  
%
% INPUT
%
% alignMovGratSt -  aligned grating structure (generated using alignMovGrtwSingleBar) 
%                   structure with .data .model and .time fields for each
%                   stim. Generated from Sandro's files using preOrganizingData
%
%
% OUTPUT
% gratingFitSt -    stracture containing moving grating data fits
%                   parameters
% gratingFitTab -   table summarizing the fits 
%
% modified after changing function that generates allCellDat to generate a
% global table. Now grating table has to be extracted


% fit parameters
numFrame = 8; % number frames in a full cycle
numCyc = 5; 
hilBuff = 100; % samples to discard from edgesof hilbert results
zeroThresh = 0.001; % cutoff for calling zero in hilbert phase result

gratingFitTab = [];
gratingFitSt = struct; 

relTab = alignMovGratSt.gratingTable;
numGStim = height(relTab);

for stim=1:numGStim

    relDur = relTab.stimDur(stim) * 1000; % to convert to ms
    
    startFudge = (0.75 + numFrame) * relDur; % to get rid of the strong onset response in slower stim 
    
    relTimeWin = [startFudge,  relDur * numFrame * (numCyc-0.5)]; % numCyc-1 to exclude the last cycle (fit just center) 

    tempTime = alignMovGratSt.result(stim).subData.baseSub(:,1); 
    tempData= alignMovGratSt.result(stim).subData.baseSub(:,2); 
    
    relInds = arrayfun(@(x) find(tempTime > x, 1, 'first'), relTimeWin); 

    relTS= tempData(relInds(1):relInds(2));
    relTime = tempTime(relInds(1):relInds(2));
    
    
    
%     plot(tempTime, tempData); hold on
    
    
    

    smoothTS = smoothdata(relTS, 'smoothingfactor', 0.1);
    meanSmTS = mean(smoothTS); 
    meanSubSmTS = smoothTS - meanSmTS; 
    
    hilRes = hilbert(meanSubSmTS); 

    resPhase = angle(hilRes);
    relResPhase = resPhase(hilBuff:end-hilBuff); 
    relTS2 = relTS(hilBuff:end-hilBuff);
    relTime2 = relTime(hilBuff:end-hilBuff); 
    minPosInd = find(abs(diff(relResPhase)) > 1.5*pi);
    maxPosInd = relResPhase > -zeroThresh & relResPhase < zeroThresh;
    
%     assert(length(minPosInd) < 6, 'too many min positions identified in stim %d - check', stim)
    
    plot(relTime2, relTS2); hold on
    line([relTime2(1), relTime2(end)], [0, 0], 'color', [1,1,1]*0.8)
    plot(relTime2, relResPhase)
    plot(relTime2(minPosInd), relTS2(minPosInd), 'o', 'color', 'r', 'markerfacecolor', 'w')
    plot(relTime2(maxPosInd), relTS2(maxPosInd), 'o', 'color', 'b', 'markerfacecolor', 'k')
    
    hold off
    
    pause
    
        
        

    
%     gratingFitSt(stim).Fit = tempFit;
%     gratingFitSt(stim).FitRes = tempRFit;
%     gratingFitSt(stim).GoF = tempGoF;
%     gratingFitSt(stim).GoFRes = tempRGoF;
% %     gratingFitSt(stim).FitOutput = tempOut;
%     gratingFitSt(stim).RelTime = relTime; 
% 
% 
%     rowTab = relTab(stim, :); 
%     meanR = tempFit.a; 
%     genSlope = tempFit.b; 
%     cosAmp = tempFit.c; 
%     cosPhase = tempFit.d; 
%     rSq = tempGoF.rsquare;
%     gratingFitTab = [gratingFitTab; [rowTab, table(meanR, genSlope, cosAmp, cosPhase, rSq, resAmp, resRSq)]]; 




end




end
