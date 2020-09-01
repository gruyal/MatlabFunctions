function varargout = plotSPFRFrom2ModelsByIter(cellNum, modelFuncName, modelDirAndIterSt, ncFlag)

% function plotSPFRFrom2ModelsByIter(cellNum, iter, ncFlag, modelString)
%
% This function is a modification of plotSPFRModelByIter. And is designed for plotting results from different fits of the same model. 
% This function plots the modeling results (from Pablo) together with the original data
% Based on plotMovModelByIter 
% 
% Inputs
%
% modelDirAndIterSt -   Struct. Contians Several fields:
% .modelDir -           String. model diretory and type where model results are saved
% .iter -               Integer. model iteration to implement from the above results 
%                      
% cellNum-              for T4 17 or 18 for T5 19 or 20 (nondiagonal and
%                       diagonal in both cases)
% iter -                1XN vector. model iterations to present (will be
%                       overlaid on top of data). one iteration per
%                       modelDir
% ncFlag -              (optional) logical. Flag to indicate if
%                       non-preferred contrast is used (defualt)
% modelFuncName -       String. model function name. can be different from
%                       model type the function.
% 
% Note!!! Use "" for modelDirAndType - otherwise split does not work. And
% no trailing /


if nargin < 4
    ncFlag = 1;
end

assert(ismember(ncFlag, [0,1]), 'ncFlag should be logical')

defMDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/modelFiles';

modelFH = str2func(modelFuncName);

allWid = [1,2,4];
allDur = [40,160];

cellT = str2double(modelFuncName(2));

if cellT == 4
    dataDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/dataFiles/T4_spfr_structs/';
    
    if ncFlag
        valSt = 'val0';
    else
        valSt = 'val1';
    end
    
    assert(ismember(cellNum, [17,18]), 'for T4 cellNum should be either 17 or 18')
    
elseif cellT == 5
    dataDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/dataFiles/T5_spfr_structs/';
    
    if ncFlag
        valSt = 'val1';
    else
        valSt = 'val0';
    end
    
    assert(ismember(cellNum, [19,20]), 'for T5 cellNum should be either 19 or 20')
end

relParams = [];
iters = zeros(length(modelDirAndIterSt), 1);
for ii=1:length(modelDirAndIterSt)

    load(fullfile(defMDir, modelDirAndIterSt(ii).modelDir), 'data_struct')
    paramTab = data_struct{cellNum}.T; 

    [~, tempParams] = extract_params2(paramTab(modelDirAndIterSt(ii).iter, :)); 
    relParams = [relParams; tempParams];
    iters(ii) = modelDirAndIterSt(ii).iter;
end



allFiles = dir(dataDir);
allFileNames = {allFiles.name};


index = 0;
spfrFilesI = zeros(length(allFileNames), 1);
spfrCondTab = table;
for ww=1:length(allWid)
    for dd=1:length(allDur)
            
            width = allWid(ww);
            duration = allDur(dd);
            
            index=index+1;
            relExp = ['spfr\w*cell', num2str(cellNum), ...
                      '_dur', num2str(duration), '_width', num2str(width), ...
                      '\w*', valSt];
                     
            tempI = (double(cellfun(@(x) ~isempty(regexp(x, relExp, 'once')), allFileNames)) * index)';
            spfrFilesI = spfrFilesI + tempI; 
            spfrCondTab = [spfrCondTab; table(index, width, duration)];
        
    end
end


prePosCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -0.03, 0.02, length(allDur) * length(allWid));
axhCell = cell(size(prePosCell));

pCol = cbrewer('qual', 'Set1', length(iters));
pCol = [[1,1,1] * 0.6; pCol];


xxLim = [0, 1000]; 

if ncFlag
    yyLim = [-7.5, 12.5];
else
    yyLim = [-10, 25];
end

figure('position', [30, 400, 1600, 700])

for ii=1:height(spfrCondTab)
    
    relFN = allFileNames{spfrFilesI == ii};
    load(fullfile(dataDir, relFN))


    numPos = length(spfr_ds.pos_vect);

    if numPos == 0
        continue
    end

    posCell = generatePositionCell(prePosCell{ii}(1), sum(prePosCell{ii}([1,3])), prePosCell{ii}(2), sum(prePosCell{ii}([2,4])), ...
            0.01, -0.1, numPos); %- 2); % I get rid of the 2 extreme positions
    axh = gobjects(size(posCell));


    % parsing position responses
    stimStInd = find(spfr_ds.stimIdx == 1);
    stimInds = [stimStInd, [stimStInd(2:end) - 1; length(spfr_ds.stimIdx)]];

    relDat = spfr_ds;

    dataVec = relDat.baseSub; 
    timeVec = relDat.time; 
    
    modelVec = zeros(length(modelDirAndIterSt), length(dataVec));
    for jj=1:length(modelDirAndIterSt)
        modelVec(jj, :) = feval(modelFH, relParams(jj,:), relDat);
    end


    for pp=1:size(stimInds, 1)

        axh(pp) = axes('position', posCell{pp});
        hold on 
        axh(pp).YLim = yyLim;
        axh(pp).XLim = xxLim;

        if pp == 1
            title(['Dur:', num2str(spfrCondTab.duration(ii)), ...
              'Wid:', num2str(spfrCondTab.width(ii))]);
        else
            axh(pp).YColor = 'none';
        end

        if ii < height(spfrCondTab)
            axh(pp).XColor = 'none';
        end

        plotDat = dataVec(stimInds(pp,1):stimInds(pp,2));
        plotMod = modelVec(:, stimInds(pp,1):stimInds(pp,2));
        plotTim = timeVec(stimInds(pp,1):stimInds(pp,2)) - timeVec(stimInds(pp,1));
        
        plot(plotTim, plotDat, 'linewidth', 1, 'color', pCol(1, :))
        for jj=1:length(iters)
            plot(plotTim, plotMod(jj,:), 'linewidth', 2, 'color', pCol(jj+1, :))
        end
        
        hold off

    end

    axhCell{ii} = axh;


end
    

preLeg = arrayfun(@num2str, iters, 'uniformoutput', 0);
legLab = cell(1, length(iters)+1);
for ll = 1:length(modelDirAndIterSt)
    legLab{ll+1} = preLeg{ll};
end
legLab{1} = 'Data';
legend(axhCell{end}(1), legLab)

    
if nargout == 1
    varargout{1} = axhCell; 
end





end