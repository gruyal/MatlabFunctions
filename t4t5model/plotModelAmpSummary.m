function varargout = plotModelAmpSummary(modSt, relIters)

% function varargout = plotModelSummary(modSt, relIters)
% This function plots a summary of several iterations from a particular
% modelData and model function
%
% INPUTS
%
% modSt -       Structure with the following fields
%   .cellType   single number (4 or 5)
%   .cellNum    single number (17 or 18 for T4; 19 or 20 for T5)
%   .modName    string. name of the model results file
%   .modFunc    string. name of model function to be used (should be in the
%               path)
%
% relIters -    vector. iteration number to be plotted


modDir  = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/';

cellT = modSt.cellType;
cellNum = modSt.cellNum;



modName = modSt.modName;

optName = ['modelFiles/', modName];

load(fullfile(modDir, optName), 'data_struct');
optTab = data_struct{cellNum}.T;


[parNames, relParams] = extract_params2(optTab(relIters, :)); 


% claculating delta Amplitude
condName = {'E_', 'I_', 'E2', 'I2'};
lineVars = {'A', 'm', 'b'};

lineVarInd = zeros(length(condName), length(lineVars));

for ii=1:length(condName)
    for jj=1:length(lineVars)
        lineVarInd(ii,jj) = find(startsWith(parNames, [lineVars{jj}, condName{ii}])); 
    end 
end

xxTime = [0, 40, 80, 160, 320, 640];

condAmpVec = nan(length(relIters), length(condName), length(xxTime));

for ii=1:length(relIters)
    for jj=1:length(condName)
        tempPar = relParams(ii, lineVarInd(jj, :)); 
        condAmpVec(ii,jj,:) = max(0, tempPar(1) * (tempPar(2) * 0.01 * xxTime + tempPar(3)));
    end
end

xSt=0.075; xEnd=0.95; ySt=0.075; yEnd=0.925;
posCell = generatePositionCell(xSt, xEnd, ySt, yEnd, 0.05, -0.015, length(condName));
axh = gobjects(size(posCell));

figure('position', [300, 700, 1000, 250])

pCol = cbrewer('seq', 'BuPu', length(relIters) + 1);
pCol = flipud(pCol);

for ii=1:length(condName)
    
    axh(ii) = axes('position', posCell{ii}, 'ColorOrder', pCol, 'NextPlot', 'replacechildren');
    
    plot(xxTime, squeeze(condAmpVec(:,ii,:)), 'linewidth', 2)
    title(condName{ii})
    axh(ii).XTick = xxTime; 
end
    




end