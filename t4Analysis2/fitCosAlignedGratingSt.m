function [gratingFitSt, gratingFitTab] = fitCosAlignedGratingSt(alignMovGratSt)

% function [gratingFitSt, gratingFitTab] = fitCosAlignedGratingSt(pStruct)
%
% This function is a modification of fitCosGratingDatAndModel using the same
% prinipels directly on the grating data (before model fitting) and meant to fit cosine to 
% the grating data and model results from Sandro. 
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
smoothWin= 25; % Sandro downsampled data
numFrame = 8; % number frames in a full cycle
numCyc = 5; 

gratingFitTab = [];
gratingFitSt = struct; 

relTab = alignMovGratSt.gratingTable;
numGStim = height(relTab);

for stim=1:numGStim

    relDur = relTab.stimDur(stim) * 1000; % to convert to ms
    
    startFudge = (0.75 + numFrame) * relDur; % to get rid of the strong onset response in slower stim 
    
    relTimeWin = [startFudge,  relDur * numFrame * (numCyc-0.5)]; % numCyc-1 to exclude the last cycle (fit just center) 
    
    stimInd = relTab.index(stim); 

    tempTime = alignMovGratSt.result(stim).subData.baseSub(:,1); 
    tempData= alignMovGratSt.result(stim).subData.baseSub(:,2); 
    
    relInds = arrayfun(@(x) find(tempTime > x, 1, 'first'), relTimeWin); 

    relTS= tempData(relInds(1):relInds(2));
    relTime = tempTime(relInds(1):relInds(2));
    
    %uncomment to see the data selected for fitting (avoiding start and
    %end)
%     plot(tempTime, tempData); hold on
%     plot(relTime, relTS); hold off
%     
%     pause

    smoothTS = smoothdata(relTS, 'movmean', smoothWin);

    stPoint = (max(smoothTS) - min(smoothTS))/2;
    fOpt = fitoptions('method', 'nonlinearleastsquare', ...
                      'startpoint', [mean(smoothTS), 0.05, stPoint, 0], ...
                      'lower', [-20, -1, 0.1, -pi], ...
                      'upper', [20, 1, 10, pi]);

    cosFit = fittype(['a + b*x + c*cos(d + 2*pi/', num2str(relDur*numFrame), '*x)'], 'options', fOpt);
    [tempFit, tempGoF, tempOut] = fit(relTime, smoothTS, cosFit);
    
    % some cases the fit doesnt capture the entire amplitude
    if relDur > 40 && tempGoF.adjrsquare < 0.9 % happens in slower grating but not all of them
        resDat = tempOut.residuals;
        
        stPoint = (max(resDat) - min(resDat))/2;
        fOpt = fitoptions('method', 'nonlinearleastsquare', ...
                          'startpoint', [mean(resDat), 0.05, stPoint, 0, relDur*(numFrame+1)], ... % added one to bias toward slower solution 
                          'lower', [-20, -1, 0.1, -pi, relDur/5], ...
                          'upper', [20, 1, 10, pi, relDur * 15]);
        resFit = fittype(['a + b*x + c*cos(d + (2*pi/f)*x)'], 'options', fOpt); % added freq as a parameter

        [tempRFit, tempRGoF] = fit(relTime, resDat, resFit);
        resAmp = tempRFit.c; 
        resRSq = tempRGoF.adjrsquare; 
    else
        
        tempRFit = [];
        tempRGoF = [];
        resAmp = NaN;
        resRSq = NaN; 
    end
        

    
    gratingFitSt(stim).Fit = tempFit;
    gratingFitSt(stim).FitRes = tempRFit;
    gratingFitSt(stim).GoF = tempGoF;
    gratingFitSt(stim).GoFRes = tempRGoF;
%     gratingFitSt(stim).FitOutput = tempOut;
    gratingFitSt(stim).RelTime = relTime; 


    rowTab = relTab(stim, :); 
    meanR = tempFit.a; 
    genSlope = tempFit.b; 
    cosAmp = tempFit.c; 
    cosPhase = tempFit.d; 
    rSq = tempGoF.rsquare;
    gratingFitTab = [gratingFitTab; [rowTab, table(meanR, genSlope, cosAmp, cosPhase, rSq, resAmp, resRSq)]]; 




end




end
