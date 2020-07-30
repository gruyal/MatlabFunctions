function modMovBarSt = calcMovModelByIter(modelDirAndType, cellNum, iter, ncFlag)

% function calcMovModelByIter(cellT, iter)
%
% This function is a modification of the plotting function generate
% the modeling results (from Pablo) together
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
%
% OUTPUT
% modMovBarSt -         structure with table and results for each condition
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
modelFH = str2func(modelType);

modMovBarSt = struct;


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
[paramN, relParams] = extract_params(paramTab(iter, :)); 


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



for ii=1:height(mvCondTab)
    
    relFN = allFileNames{movFilesI == ii};
    load(fullfile(dataDir, relFN))
    
    if mvCondTab.direction(ii) == 1
        relDat = pd_ds;
    else
        relDat = nd_ds;
    end
    
    modelVec = feval(modelFH, relParams, relDat);
    dataVec = relDat.baseSub; 
    timeVec = relDat.time; 
    
    modMovBarSt.result(ii).model = modelVec;
    modMovBarSt.result(ii).data = dataVec;
    modMovBarSt.result(ii).time = timeVec;
    
    
end


modMovBarSt.table = mvCondTab;


end