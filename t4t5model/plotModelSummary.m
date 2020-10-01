
function varargout = plotModelSummary(modSt, relIters)

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

if cellT == 4
    dataDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/dataFiles/T4_spfr_structs/';
    assert(ismember(cellNum, [17,18]), 'for T4 cellNum should be either 17 or 18')
    ncVal = 0;
    pcVal = 1;
elseif cellT == 5
    dataDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/dataFiles/T5_spfr_structs/';
    assert(ismember(cellNum, [19,20]), 'for T5 cellNum should be either 19 or 20')
    ncVal = 1;
    pcVal = 0;
end

modName = modSt.modName;
modFunc = modSt.modFunc;
modelFH = str2func(modFunc);

optName = ['modelFiles/', modName];

allFiles = dir(fullfile(dataDir));
allFileNames = {allFiles.name};

load(fullfile(modDir, optName), 'data_struct');
optTab = data_struct{cellNum}.T;

varNames = optTab.Properties.VariableNames;
[parNames, relParams] = extract_params2(optTab(relIters, :)); 
tausInds = startsWith(varNames, 'T');
relTaus = optTab{relIters, tausInds} * 10; % organized E, I, E2, I2 (rise, decay) and upstream
upTauNI = startsWith(parNames, 'Ti');
upTauN = cellfun(@(x) x(1:4), parNames(upTauNI), 'UniformOutput', false);

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


allWid = 4;
allDur = [40,160];
allDir = {'p', 'n'};
allVal = [0,1];

index = 0;
movFilesI = zeros(length(allFileNames), 1);
mvCondTab = table;
for ww=1:length(allWid)
    for dd=1:length(allDur)
        for dr=1:length(allDir)
            for vv=1:length(allVal)
            
                if eq(allDir{dr}, 'p')
                    direction = 1;
                else
                    direction = 0;
                end

                width = allWid(ww);
                duration = allDur(dd);
                val = allVal(vv);

                index=index+1;
                relExp = [allDir{dr},'d\w*cell', num2str(cellNum), ...
                          '_dur', num2str(duration), '_width', num2str(width), ...
                          '\w*', 'val', num2str(val)];

                tempI = (double(cellfun(@(x) ~isempty(regexp(x, relExp, 'once')), allFileNames)) * index)';
                movFilesI = movFilesI + tempI; 
                mvCondTab = [mvCondTab; table(index, width, duration, direction, val)];
                
            end
        end
    end
end

% for the last row 
index = 0;
movFilesI2 = zeros(length(allFileNames), 1);
mvCondTab2 = table;
for dd=1:length(allDur)
    for dr=1:length(allDir)
        

        if eq(allDir{dr}, 'p')
            direction = 1;
        else
            direction = 0;
        end

        width = 2;
        duration = allDur(dd);
        val = ncVal;

        index=index+1;
        relExp = [allDir{dr},'d\w*cell', num2str(cellNum), ...
                  '_dur', num2str(duration), '_width', num2str(width), ...
                  '\w*', 'val', num2str(val)];

        tempI = (double(cellfun(@(x) ~isempty(regexp(x, relExp, 'once')), allFileNames)) * index)';
        movFilesI2 = movFilesI2 + tempI; 
        mvCondTab2 = [mvCondTab2; table(index, width, duration, direction, val)];

        
    end
end



condMax = zeros(size(relParams,1), height(mvCondTab), 6); % 4 conductances + 2 max total
muNSig = zeros(size(relParams,1), 2, 4); % 4 condunctances, mu and sigma
maxDNM = zeros(size(relParams,1), height(mvCondTab), 2);
vTraces = cell(size(relParams,1), height(mvCondTab2), 2);
timeVecs = cell(height(mvCondTab2), 1);

for ii=1:size(relParams,1)
    for jj=1:height(mvCondTab)
    
        relFN = allFileNames{movFilesI == jj};
        load(fullfile(dataDir, relFN))

        if mvCondTab.direction(jj) == 1
            relDat = pd_ds;
        else
            relDat = nd_ds;
        end

        [vVec, totGE, totGI, gMat, sMat] = feval(modelFH, relParams(ii,:), relDat);
        condMax(ii,jj, :) = [cellfun(@max, gMat), max(totGE), max(totGI)];
        maxDNM(ii, jj, 1) = quantile(relDat.baseSub, 0.995);
        maxDNM(ii, jj, 2) = quantile(vVec, 0.995);
    end
    
    muNSig(ii, :, :) = reshape([sMat{:}], 2, []); 
    
    for jj=1:height(mvCondTab2)
    
        relFN2 = allFileNames{movFilesI2 == jj};
        load(fullfile(dataDir, relFN2))

        if mvCondTab2.direction(jj) == 1
            relDat = pd_ds;
        else
            relDat = nd_ds;
        end

        vVec = feval(modelFH, relParams(ii,:), relDat);
        
        vTraces{ii, jj, 1} = relDat.baseSub;
        vTraces{ii, jj, 2} = vVec;
        timeVecs{jj} = relDat.time;
    end
    
    
end


relCol = cbrewer('qual', 'Paired', 10);
relCol = relCol([6,2,5,1], :); % to match the order of gMat 
modCol = [[1,1,1]*0.75; [1,1,1]*0.8];

xSt=0.075; xEnd=0.95; ySt=0.05; yEnd=0.875;
posCellU = generatePositionCell(xSt, xEnd, yEnd+0.025, yEnd+0.1 , 0.01, -0.075, length(relIters));
posCell = generatePositionCell(xSt, xEnd, ySt, yEnd, 0.01, 0.015, [length(relIters), 6]);
axh = gobjects(size(posCell) + [0,1]);

figure('position', [40, 60, 1500, 900])

condLab = {'E', 'I', 'E2', 'I2'};
tausLab = [{'Tr'}, {'Td'}, upTauN(:)'];
stimLab1 = {'Fast', 'Fast', 'Slow', 'Slow'}; % based on order in mvCondTab
stimLab2 = {'PD', 'ND', 'PD', 'ND'}; 
stimLabArray = [stimLab1; stimLab2];
stimTickLab = strtrim(sprintf('%s\\newline%s\n', stimLabArray{:}));
edgeCol = [0,0,0; [1,1,1]*0.8];

for ii=1:length(relIters)
    
    axh(ii,1) = axes('position', posCellU{ii,1}, 'xlim', [-6.5, 6.5], 'ylim', [0.75, 4.25]); 
    hold on 
    for jj=1:size(muNSig,3)
        plot(muNSig(ii,1,jj), jj, 'o', 'markeredgecolor', relCol(jj,:), ...
             'markerfacecolor', relCol(jj,:), 'markersize', 8)
        line([muNSig(ii,1,jj) - muNSig(ii,2,jj), muNSig(ii,1,jj) + muNSig(ii,2,jj)], ones(1,2) * jj, ...
             'linewidth', 4, 'color', relCol(jj,:))      
    end
    hold off
    axh(ii,1).YTick = 1:length(condLab);
    if ii==1
        axh(ii,1).YTickLabel = condLab;
    else
        axh(ii,1).YTickLabel = {};
    end
    title(num2str(relIters(ii)))
    
    
    
    axh(ii,2) = axes('position', posCell{ii,1}, 'xlim', [0.75, size(relTaus,2)/3+0.5]);
    hold on 
    
    for jj=1:size(muNSig,3)
        plot([1,2], [relTaus(ii,jj), relTaus(ii,jj+4)] , 'linewidth', 4, ...
             'color', relCol(jj,:), 'marker', 'o', 'markeredgecolor', relCol(jj,:), 'markerfacecolor', 'w')
    end
    
    for jj=9:size(relTaus, 2) % since simple and complex have different upstream filters
        plot(jj/3, relTaus(ii,jj) , ...
             'color', 'k', 'marker', 'o', 'markeredgecolor', 'k', 'markerfacecolor', 'k')
    end
    hold off
    axh(ii,2).XTick = [1,2, (9:size(relTaus, 2))/3];
    axh(ii,2).XTickLabel = tausLab;
    xtickangle(45)
    
    if ii>1
        axh(ii,2).YTickLabel = {};
    end
    
    axh(ii,3) = axes('position', posCell{ii,2}, 'xlim', [0.5, length(stimLab1)+0.5]);
    hold on
    pcInd = mvCondTab.val == pcVal; 
    for jj=1:size(muNSig,3)
        plot(1:sum(pcInd), squeeze(condMax(ii,pcInd, jj)), 'linewidth', 4, 'color', relCol(jj,:), ...
             'marker', 'o', 'markeredgecolor', relCol(jj,:), 'markerfacecolor', relCol(jj,:))
    end
    axh(ii,3).XTick = 1:length(stimLab1);
    axh(ii,3).XTickLabel = {};
    
    if ii==1
        axh(ii,3).YLabel.String = {'Max Curr by C'; 'PC Stim'};
    else
        axh(ii,3).YTickLabel = {};
    end
    
    hold off
    
    axh(ii,4) = axes('position', posCell{ii,3}, 'xlim', [0.5, length(stimLab1)+0.5]);
    hold on
    ncInd = mvCondTab.val == ncVal; 
    for jj=1:size(muNSig,3)
        plot(1:sum(ncInd), squeeze(condMax(ii,ncInd, jj)), 'linewidth', 4, 'color', relCol(jj,:), ...
             'marker', 'o', 'markeredgecolor', relCol(jj,:), 'markerfacecolor', relCol(jj,:))
    end
    hold off
    
    axh(ii,4).XTick = 1:length(stimLab1);
    axh(ii,4).XTickLabel = stimTickLab;
    
    if ii==1
        axh(ii,4).YLabel.String = {'Max Curr by C'; 'NC Stim'};
    else
        axh(ii,4).YTickLabel = {};
    end
    
    
    axh(ii,5) = axes('position', posCell{ii,4}, 'xlim', [0.5, size(maxDNM,2)+0.5]);
    hold on
    
    for jj=1:size(maxDNM,2)
        
        markEC = 1 + mod(jj-1,2);
        markS = 4 + (jj > 4) * 4;
        markFC = 1+ismember(jj,[3,4,7,8]);
        
        plot(ones(2,1) * jj, squeeze(maxDNM(ii,jj,:)), 'linewidth', 1, 'color', 'k')
        plot(jj, maxDNM(ii,jj,1),  'marker', 'o', 'markeredgecolor', edgeCol(markEC,:), 'markerfacecolor', edgeCol(markEC,:))
        plot(jj, maxDNM(ii,jj,2),  'marker', 'o', 'markeredgecolor', 'k', 'markerfacecolor', relCol(markFC,:), 'markersize', markS)
             
    end
    hold off
    
    axh(ii,5).XTick = 1:size(maxDNM,2);
    axh(ii,5).XTickLabel = {};
    
    if ii==1
        axh(ii,5).YLabel.String = {'Max V'; 'Data v. mod'};
    else
        axh(ii,5).YTickLabel = {};
    end
    
    axh(ii,6) = axes('position', posCell{ii,5}, 'xlim', [0, xxTime(end)]);
    hold on
    
    for jj=1:size(condAmpVec,2)
        line([0, xxTime(end)], [0,0], 'color', 'k', 'linewidth', 1)
        plot(xxTime, squeeze(condAmpVec(ii,jj,:)), 'linewidth', 2, 'color', relCol(jj,:))
        
    end
    hold off
    
    axh(ii,6).XTick = [40, 160, 320, 640]; 
    axh(ii,6).XTickLabel = arrayfun(@num2str,[40, 160, 320, 640], 'uniformoutput', 0);
    
    if ii==1
        axh(ii,6).YLabel.String = {'amp v. effDur'};
    end
    
    axh(ii,7) = axes('position', posCell{ii,6});%, 'xlim', [0.5, size(maxDNM,2)+0.5]);
    hold on
    
    for jj=1:size(vTraces,2)
        
        markC = 1 + mod(jj-1,2);
        addV = (jj < 3) * 10; 
        
        plot(timeVecs{jj}, vTraces{ii,jj,1} + addV, 'linewidth', 1, 'color', modCol(markC,:))
        plot(timeVecs{jj}, vTraces{ii,jj,2} + addV, 'linewidth', 2, 'color', relCol(markC,:))
             
    end
    hold off
    
    
    axh(ii,7).XTickLabel = {};
    
    if ii==1
        axh(ii,7).YLabel.String = {'Response W2 NC'; 'Data v. mod'};
    else
        axh(ii,7).YTickLabel = {};
    end
    
    
end


axSz = size(axh);
if axSz(1) > 1
    for jj=2:axSz(2)
        if jj==6 % since scales can be very different
            continue
        end
        tempRg = zeros(axSz(1),2); 
        for ii=1:axSz(1)
            tempRg(ii, :) = axh(ii,jj).YLim;
        end
        maxRg = max(tempRg);
        for ii=1:axSz(1)
            axh(ii,jj).YLim = maxRg;
        end
    end
end

if nargout ==1
    varargout{1} = axh;
end



end