function [gratingFitSt, gratingFitTab] = fitCosGratingDatAndModel(allCellDat, fitAllDatFlag)

% function gratingFitSt = fitCosGratingDatAndModel(pStruct)
%
% This function is a modification of fitCosFlickerProtocolPiecewise and meant to fit cosine to 
% the grating data and model results from Sandro. 
% As a poor man version of de-trending, it performs peicewise fits and
% averages the response
%
% % change function since mgTab is now included in the allCellDat structure
%
% INPUT
%
% allCellDat -      structure with .data .model and .time fields for each
%                   stim. Generated from Sandro's files using preOrganizingData
% fitAllDatFlag -   (optional) logical. If true cos fit is done on all stim
%                   if false, fit only the long duration stimuli. defualt
%                   is true
%
%
% OUTPUT
% gratingFitSt -    stracture containing moving grating data and model fits
% gratingFitTab -   table summarizing the fits 
%
% modified after changing function that generates allCellDat to generate a
% global table. Now grating table has to be extracted

if nargin < 2 
    fitAllDatFlag = 1; 
end

longDur = 160; 

% fit parameters
smoothWin= 25; % Sandro downsampled data
numFrame = 8; % number frames in a full cycle

gratingFitTab = [];


if fitAllDatFlag
    relTab = allCellDat.table(allCellDat.table.protType == 1, :);
    numGStim = height(relTab);
else
    relTab = allCellDat.table(allCellDat.table.protType == 1 & allCellDat.table.duration == longDur, :);
    numGStim = height(relTab);
end

for stim=1:numGStim

    relDur = relTab.duration(stim); 
    if relDur > 40
        startFudge = 500; % to get rid of the strong onset response in slower stim 
    else
        startFudge = 100; 
    end
    relTimeWin = [startFudge,  relDur * 25 - startFudge]; % 25 is a estimate of when stim ends
    
    stimInd = relTab.index(stim); 

    tempTime = allCellDat.MG(stimInd).time; 
    tempData= allCellDat.MG(stimInd).data; 
    tempModel= allCellDat.MG(stimInd).model;         
%         [cc, stim]
    relInds = arrayfun(@(x) find(tempTime > x, 1, 'first'), relTimeWin); 

    relTS= tempData(relInds(1):relInds(2));
    relModTS= tempModel(relInds(1):relInds(2));
    relTime = tempTime(relInds(1):relInds(2));

    tempD = {relTS, relModTS};

    for dd = 1:length(tempD)

        smoothTS = smooth(tempD{dd}, smoothWin);

        stPoint = (max(smoothTS) - min(smoothTS))/2;
        fOpt = fitoptions('method', 'nonlinearleastsquare', ...
                          'startpoint', [mean(smoothTS), 0.05, stPoint, 0], ...
                          'lower', [-20, -1, 0.1, -pi], ...
                          'upper', [20, 1, 10, pi]);

        cosFit = fittype(['a + b*x + c*cos(d + 2*pi/', num2str(relDur*numFrame), '*x)'], 'options', fOpt);

        [tempFit, tempGoF] = fit(relTime, smoothTS, cosFit);

        if dd==1
            gratingFitSt.MG(stim).dataFit = tempFit;
            gratingFitSt.MG(stim).dataGof = tempGoF;
            gratingFitSt.MG(stim).timeFit = relTime; 
            dataFlag = 1; 
        else
            gratingFitSt.MG(stim).modelFit = tempFit;
            gratingFitSt.MG(stim).modelGof = tempGoF;
            dataFlag = 0; 
        end


        rowTab = relTab(stim, :); 
        meanR = tempFit.a; 
        genSlope = tempFit.b; 
        cosAmp = tempFit.c; 
        cosPhase = tempFit.d; 
        rSq = tempGoF.rsquare;
        gratingFitTab = [gratingFitTab; [rowTab, table(meanR, genSlope, cosAmp, cosPhase, rSq, dataFlag)]]; 


    end

end




end
