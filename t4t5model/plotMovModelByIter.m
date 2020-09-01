function varargout = plotMovModelByIter(modelDirAndType, cellNum, iter, ncFlag, modelString)

% function plotMovModelByIter(cellT, iter)
%
% This function plots the modeling results (from Pablo) together with the original data
%
% Inputs
%
% modelDirAndType -     String. model diretory and type (file name) which includes the the cell type (T4 or
%                       T5) and the type of model implemented
%                       (classic/delta/etc). Directory is appended to the
%                       defDir
% cellNum-              for T4 17 or 18 for T5 19 or 20 (nondiagonal and
%                       diagonal in both cases)
% iter -                integer. model iteration to present
% ncFlag -              (optional) logical. Flag to indicate if
%                       non-preferred contrast is used (defualt)
% modelString -         (optional) if model function name is different from
%                       model type the function can be entered here as a string
% 
% 
% Note!!! Use "" for modelDirAndType - otherwise split does not work. And
% no trailing /


if nargin < 4
    ncFlag = 1; 
end

assert(ismember(ncFlag, [0,1]), 'ncFlag should be logical')

defMDir = '/Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/T4recordingSummaryAndAnalysis2/pabloModAna/modelFiles';

load(fullfile(defMDir, modelDirAndType), 'data_struct')
spRes = split(string(modelDirAndType), "/");

modelType = spRes{end};
if nargin < 5
    modelFH = str2func(modelType);
else
    modelFH = str2func(modelString);
end


allWid = [1,2,4];
allDur = [40,160];
allDir = {'p', 'n'};

cellT = str2double(modelType(2));

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

paramTab = data_struct{cellNum}.T; 
[paramN, relParams] = extract_params2(paramTab(iter, :)); 


allFiles = dir(dataDir);
allFileNames = {allFiles.name};


index = 0;
movFilesI = zeros(length(allFileNames), 1);
mvCondTab = table;
for ww=1:length(allWid)
    for dd=1:length(allDur)
        for dr=1:length(allDir)
            
            if eq(allDir{dr}, 'p')
                direction = 1;
            else
                direction = 0;
            end
            
            width = allWid(ww);
            duration = allDur(dd);
            
            index=index+1;
            relExp = [allDir{dr},'d\w*cell', num2str(cellNum), ...
                      '_dur', num2str(duration), '_width', num2str(width), ...
                      '\w*', valSt];
                     
            tempI = (double(cellfun(@(x) ~isempty(regexp(x, relExp, 'once')), allFileNames)) * index)';
            movFilesI = movFilesI + tempI; 
            mvCondTab = [mvCondTab; table(index, width, duration, direction)];
        end
    end
end


posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.03, 0.02, [length(allDur), length(allWid)]);
axh = gobjects(size(posCell)); 

pCol = cbrewer('qual', 'Paired', 6);
pCol = pCol([5,6,1,2], :);


xxLimM = [1500; 3200]; 

if ncFlag
    yyLim = [-7.5, 17.5];
else
    yyLim = [-10, 30];
end

figure('position', [750, 200, 700, 750], 'name', num2str(iter))


for ii=1:height(mvCondTab)
    
    relFN = allFileNames{movFilesI == ii};
    load(fullfile(dataDir, relFN))
    
    if mvCondTab.direction(ii) == 1
        relDat = pd_ds;
        colInd = 2; 
        axInd = ceil(ii/2);
        axh(axInd) = axes('position', posCell{axInd});
        axh(axInd).YLim = yyLim;
        
        if mvCondTab.duration(ii) == allDur(1)
            axh(axInd).XLim = [0, xxLimM(1)];
        else
            axh(axInd).XLim = [0, xxLimM(2)];
        end
        
        if axInd < 3
            title(num2str(mvCondTab.duration(ii)))
        end
        
        if mod(axInd,2) == 1
            ylabel(mvCondTab.width(ii), 'fontweight', 'bold')
        else
            axh(axInd).YColor = 'none';
        end
           
        if axInd < height(mvCondTab)/2 - length(allDur) + 1
            axh(axInd).XColor = 'none';
        end
        
        hold on 
    else
        relDat = nd_ds;
        colInd = 4; 
    end
    
    modelVec = feval(modelFH, relParams, relDat);
    dataVec = relDat.baseSub; 
    timeVec = relDat.time; 
    
    
    
    plot(timeVec, dataVec, 'linewidth', 1, 'color', pCol(colInd, :))
    plot(timeVec, modelVec, 'linewidth', 2, 'color', pCol(colInd-1, :))
    
    if colInd == 4 % nd
        hold off
    end
    
    
end


    
if nargout == 1
    varargout{1} = axh; 
end





end