function allCellDat = organizingClusterData(cellNum, iterNum, relDir)

% function allCellDat = organizingClusterData(relDir, cellNum, iterNum)
%
% This function loads all the data from one cell and one iteration of the
% model and organizes it in a structure and tables
%
% INPUT
% relDir -          (optional) relevant directory. if not given current
%                   directory is used.
% cellNum -         cell number to load
% interNum -        iteration number from the model (implemented as a
%                   directory name for that iteration)
%
% OUTPUT
% allCellDat -      structure containing original data and model simulation
%                   and relevant tables that reference them 


%{
0: single bars
1: moving grating
2: static grating
3: moving bars
4: minimal motion 
%} 

vCutoff = 2.5; % cutoff in mV under which error is not calculated (applied to all stimuli evenly) 

if nargin < 3
    relDir = pwd; 
end

protOrd = {'SB'; 'MG'; 'SG'; 'MB'; 'MM'};

sbTable = [];
mbTable = [];
sgTable = [];
mgTable = [];
mmTable = [];

timeBuff = 30; % for minMot data



dataSName = fullfile(relDir, ['/data_cell_',num2str(cellNum),'_all.mat']);
load(dataSName,'p','d')    

modelSName = fullfile(relDir, num2str(iterNum),['/result_cell_',num2str(cellNum),'.mat']);
textVName = ['cond_viol_',num2str(cellNum),'.txt'];

dirRes = dir(fullfile(relDir, num2str(iterNum)));

load(modelSName, 'vm_all','q', 'err_spfr');
d.vm = vm_all; 

allCellDat.errorSPFR = err_spfr; 
allCellDat.params = q.param;
allCellDat.fitErr = q.resnorm; 

if any(cellfun(@(x) strcmp(x, textVName), {dirRes.name}))
    violTab = readtable(fullfile(relDir, num2str(iterNum), textVName)); 
    if height(violTab) > 0 % table isnt empty
        violTab.Properties.VariableNames = {'protType', 'posOrDir', 'width', 'dur', 'prop', 'maxVal'};
        maxNormProp = (violTab.maxVal - q.param(end-1))/ q.param(end-1); % dividing by e2star
        violTab = addvars(violTab, maxNormProp);
    end
else
    violTab = [];
end


relProts = unique(p.protocol);

indOffset = 0; 

% single bar 

if ismember(0, relProts)
    
    sbNumStim = length(d.ind_ref1_sb);
    timeInd = p.protocol == 0; 
    relTime = p.t(timeInd); 
    maxData = nan(sbNumStim,1);
    maxMod = nan(sbNumStim,1);
    dataLen = nan(sbNumStim,1);
    fracErrSS = nan(sbNumStim,1);
    fullErrSS = nan(sbNumStim,1);
    protInd = find(cellfun(@(x) strcmp('SB',x),  protOrd));
    protType = ones(sbNumStim,1) * (protInd-1);
    direction  = maxMod; phase = maxMod; 
    
    for ii=1:sbNumStim
        relInds = d.ind_ref1_sb(ii):d.ind_ref2_sb(ii); 
        allCellDat.SB(ii).data = d.vm_all_sb(relInds); 
        allCellDat.SB(ii).model= d.vm(relInds + indOffset); 
        allCellDat.SB(ii).time = relTime{ii};
        tempData = d.vm_all_sb(relInds);
        tempMod = d.vm(relInds + indOffset);
        maxData(ii) = max(tempData); 
        maxMod(ii) = max(tempMod);
        tempErr = ((tempData - tempMod)./tempData).^2; 
        tempFullEr = (tempData - tempMod).^2; 
        relErrInd = tempData > vCutoff; 
        fracErrSS(ii) = sum(tempErr(relErrInd)); 
        if sum(relErrInd) %sorrecting for the length of the calculation 
            fullErrSS(ii) = sum(tempFullEr(relErrInd)) / sum(relErrInd);
        else
            fullErrSS(ii) = NaN;
        end
        dataLen(ii) = length(relInds); 
    end

    tempTab = table(protType, (1:sbNumStim)', p.location_sb', p.duration_sb', p.width_sb', maxData, maxMod, fracErrSS, fullErrSS, dataLen, direction, phase);
    tempTab.Properties.VariableNames = {'protType'; 'index'; 'position'; 'duration'; 'width'; 'maxData'; 'maxMod'; 'fracError'; 'fullError'; 'dataLength'; 'direction'; 'phase'};
    sbTable = vertcat(sbTable, tempTab); 

    clear tempTab maxData maxMod dataLen

    indOffset = indOffset + d.ind_ref2_sb(end); 

end

% moving bar

if ismember(3, relProts)

    mvNumStim = length(d.ind_ref1_mb);
    timeInd = p.protocol == 3; 
    relTime = p.t(timeInd); 
    maxData = nan(mvNumStim,1);
    maxMod = nan(mvNumStim,1);
    dataLen = nan(mvNumStim,1);
    fracErrSS = nan(mvNumStim,1);
    fullErrSS = nan(mvNumStim,1);
    protInd = find(cellfun(@(x) strcmp('MB',x),  protOrd));
    protType = ones(mvNumStim,1) * (protInd-1);
    phase = maxMod; position = maxMod; 

    for ii=1:mvNumStim
        relInds = d.ind_ref1_mb(ii):d.ind_ref2_mb(ii); 
        allCellDat.MB(ii).data = d.vm_all_mb(relInds); 
        allCellDat.MB(ii).model= d.vm(relInds + indOffset); 
        allCellDat.MB(ii).time = relTime{ii};
        tempData = d.vm_all_mb(relInds); 
        tempMod = d.vm(relInds + indOffset);
        maxData(ii) = max(tempData);
        maxMod(ii) = max(tempMod);
        tempErr = ((tempData - tempMod) ./ tempData).^2; 
        tempFullEr = (tempData - tempMod).^2; 
        relErrInd = tempData > vCutoff; 
        fracErrSS(ii) = sum(tempErr(relErrInd)); 
        if sum(relErrInd) %sorrecting for the length of the calculation 
            fullErrSS(ii) = sum(tempFullEr(relErrInd)) / sum(relErrInd);
        else
            fullErrSS(ii) = NaN;
        end
        dataLen(ii) = length(relInds); 
    end

    tempTab = table(protType, (1:mvNumStim)', p.direction_mb', p.duration_mb', p.width_mb', maxData, maxMod, fracErrSS, fullErrSS, phase, position, dataLen);
    tempTab.Properties.VariableNames = {'protType'; 'index'; 'direction'; 'duration'; 'width'; 'maxData'; 'maxMod'; 'fracError'; 'fullError'; 'phase'; 'position'; 'dataLength'};
    mbTable = vertcat(mbTable, tempTab); 

    clear tempTab maxData maxMod

    indOffset = indOffset + d.ind_ref2_mb(end); 

end

% moving grating

if ismember(2, relProts)

    mgNumStim = length(d.ind_ref1_mg);
    timeInd = p.protocol == 1; 
    relTime = p.t(timeInd); 
    maxData = nan(mgNumStim,1);
    maxMod = nan(mgNumStim,1);
    dataLen = nan(mgNumStim,1);
    fracErrSS = nan(mgNumStim,1);
    fullErrSS = nan(mgNumStim,1);
    protInd = find(cellfun(@(x) strcmp('MG',x),  protOrd));
    protType = ones(mgNumStim,1) * (protInd-1);
    position = maxData; width = maxData; 
    
    % left maxData and maxMod empty since it is a bit meanigless in this
    % stimulus
    for ii=1:mgNumStim
        relInds = d.ind_ref1_mg(ii):d.ind_ref2_mg(ii); 
        allCellDat.MG(ii).data = d.vm_all_mg(relInds); 
        allCellDat.MG(ii).model= d.vm(relInds + indOffset); 
        allCellDat.MG(ii).time = relTime{ii};
        tempData = d.vm_all_mg(relInds); 
        tempMod = d.vm(relInds + indOffset);
        tempErr = ((tempData - tempMod) ./ tempData).^2; 
        tempFullEr = (tempData - tempMod).^2; 
        relErrInd = tempData > (vCutoff / 2); % divide by 2 since grating is only stim that oscillates
        fracErrSS(ii) = sum(tempErr(relErrInd)); 
        if sum(relErrInd) %sorrecting for the length of the calculation 
            fullErrSS(ii) = sum(tempFullEr(relErrInd)) / sum(relErrInd);
        else
            fullErrSS(ii) = NaN;
        end
        dataLen(ii) = length(relInds); 
    end
    
    preDirec = diff(p.phase_mg(:,2:3), [], 2); 
    dirFlipFlag = sign(p.posD_sg{1}(end));
    direction = fix(preDirec * dirFlipFlag + 0.01); 
    if dirFlipFlag == 1
        tempPhase = p.phase_mg(:,1); 
    else
        tempPhase = flipud(p.phase_mg(:,1)); 
    end
    
    tempTab = table(protType, (1:mgNumStim)', direction, p.stimdur_mg, tempPhase, maxData, maxMod, fracErrSS, fullErrSS, dataLen, position, width);
    tempTab.Properties.VariableNames = {'protType'; 'index'; 'direction'; 'duration'; 'phase'; 'maxData'; 'maxMod'; 'fracError'; 'fullError';'dataLength'; 'position'; 'width'};
    mgTable = vertcat(mgTable, tempTab); 
    allCellDat.dir = dirFlipFlag;

    clear tempTab

    indOffset = indOffset + d.ind_ref2_mg(end); 

end


% static grating 

if ismember(1, relProts)

    sgNumStim = length(d.ind_ref1_sg);
    timeInd = p.protocol == 2; 
    relTime = p.t(timeInd); 
    maxData = nan(sgNumStim,1);
    maxMod = nan(sgNumStim,1);
    dataLen = nan(sgNumStim,1);
    fracErrSS = nan(sgNumStim,1);
    fullErrSS = nan(sgNumStim,1);
    protInd = find(cellfun(@(x) strcmp('SG',x),  protOrd));
    protType = ones(sgNumStim,1) * (protInd-1);
    position = maxData; width = maxData; direction = maxData; 
    

    for ii=1:sgNumStim
        relInds = d.ind_ref1_sg(ii):d.ind_ref2_sg(ii); 
        allCellDat.SG(ii).data = d.vm_all_sg(relInds); 
        allCellDat.SG(ii).model= d.vm(relInds + indOffset); 
        allCellDat.SG(ii).time = relTime{ii};
        tempData = d.vm_all_sg(relInds); 
        tempMod = d.vm(relInds + indOffset);
        maxData(ii) = max(tempData); 
        maxMod(ii) = max(tempMod); 
        tempErr = ((tempData - tempMod) ./ tempData).^2; 
        tempFullEr = (tempData - tempMod).^2; 
        relErrInd = tempData > vCutoff; 
        fracErrSS(ii) = sum(tempErr(relErrInd)); 
        if sum(relErrInd) %sorrecting for the length of the calculation 
            fullErrSS(ii) = sum(tempFullEr(relErrInd)) / sum(relErrInd);
        else
            fullErrSS(ii) = NaN;
        end
        dataLen(ii) = length(relInds); 
    end

    phaseInd = nan(length(p.posD_sg), 1); 
    for pp=1:length(p.posD_sg)
        tmp = cellfun(@(x) isequal(x, p.posD_sg{pp}), p.PosD_mg);
        phaseInd(pp) = find(tmp); 
    end

    tempTab = table(protType, (1:sgNumStim)', p.stimdur_sg, phaseInd, maxData, maxMod, fracErrSS, fullErrSS, dataLen, position, width, direction);
    tempTab.Properties.VariableNames = {'protType'; 'index'; 'duration'; 'phase'; 'maxData'; 'maxMod'; 'fracError'; 'fullError'; 'dataLength'; 'position'; 'width'; 'direction'};
    sgTable = vertcat(sgTable, tempTab); 

    clear tempTab maxData maxMod

    indOffset = indOffset + d.ind_ref2_sg(end);

end

% minmmal motion

if ismember(4, relProts)

    mmNumStim = length(d.ind_ref1_mm);
    timeInd = p.protocol == 4; 
    relTime = p.t(timeInd); 
    maxData = nan(mmNumStim,2); % for both first and second bars
    maxMod = nan(mmNumStim,2);
    fbTime = nan(mmNumStim,2);
    sbTime = nan(mmNumStim,2);
    dataLen = nan(mmNumStim,1);
    fracErrSS = nan(mmNumStim,1);
    fullErrSS = nan(mmNumStim,1);
    protInd = find(cellfun(@(x) strcmp('MM',x),  protOrd));
    protType = ones(mmNumStim,1) * (protInd-1);
    position = dataLen; phase = dataLen; 
    
    % change the extracted details so that they will fit in a global table
    for ii=1:mmNumStim
        relInds = d.ind_ref1_mm(ii):d.ind_ref2_mm(ii); 
        tempMMD = d.vm_all_mm(relInds); 
        tempMMM = d.vm(relInds + indOffset);
        tempT = relTime{ii};

        allCellDat.MM(ii).data = tempMMD;
        allCellDat.MM(ii).model= tempMMM;
        allCellDat.MM(ii).time = tempT; 
        fbTime(ii, :) = [p.fbton_mm, p.fbtoff_mm(ii)] +timeBuff;
        sbTime(ii, :) = [p.sbton_mm(ii), p.sbtoff_mm(ii)] + timeBuff;

        fbInd = arrayfun(@(x) find(tempT - x > 0, 1, 'first'), fbTime(ii, :));
        sbInd = arrayfun(@(x) find(tempT - x > 0, 1, 'first'), sbTime(ii, :));

        maxData(ii, 1) = max(tempMMD(fbInd(1):fbInd(2)));
        maxData(ii, 2) = max(tempMMD(sbInd(1):sbInd(2)));

        maxMod(ii, 1) = max(tempMMM(fbInd(1):fbInd(2)));
        maxMod(ii, 2) = max(tempMMM(sbInd(1):sbInd(2)));
        
        tempErr = ((tempMMD - tempMMM) ./ tempMMD).^2; 
        tempFullEr = (tempMMD - tempMMM).^2; 
        relErrInd = tempMMD > vCutoff; 
        fracErrSS(ii) = sum(tempErr(relErrInd)); 
        if sum(relErrInd) %sorrecting for the length of the calculation 
            fullErrSS(ii) = sum(tempFullEr(relErrInd)) / sum(relErrInd);
        else
            fullErrSS(ii) = NaN;
        end

        dataLen(ii) = length(relInds);

    end

%     tempTab = table(ones(mmNumStim, 1) * cellNum, (1:mmNumStim)', p.duration_mm, p.timeDiff_mm , p.width_mm, p.fb_pos_mm, p.sb_pos_mm, p.speedCor_mm, maxData, maxMod, fbTime, sbTime);
%     tempTab.Properties.VariableNames = {'cellNum'; 'index'; 'duration'; 'timeDiff'; 'width'; 'FBPos'; 'SBPos'; 'speedCor'; 'maxData'; 'maxMod'; 'fbTime'; 'sbTime'};
    preDirec = sign(p.sb_pos_mm - p.fb_pos_mm);
    preDirec(preDirec == 0) = NaN; 
    direction = fix(preDirec + 0.01); % to turn the difference into 0 ND and 1 PD vector
    tempTab = table(protType, (1:mmNumStim)', p.duration_mm , p.width_mm, direction, max(maxData, [], 2), max(maxMod, [], 2), fracErrSS, fullErrSS, dataLen, position, phase);
    tempTab.Properties.VariableNames = {'protType'; 'index'; 'duration'; 'width'; 'direction'; 'maxData'; 'maxMod'; 'fracError'; 'fullError'; 'dataLength'; 'position'; 'phase'};
    mmTable = vertcat(mmTable, tempTab); 
    
    cell4mmTab = ones(height(mmTable), 1) * cellNum;
    iter4mmTab = ones(height(mmTable), 1) * iterNum;
    
    fullMMTable = table(cell4mmTab, iter4mmTab, (1:mmNumStim)', p.fb_pos_mm, p.sb_pos_mm, p.speedCor_mm, p.timeDiff_mm, p.width_mm, p.duration_mm, maxData(:,1), maxMod(:,1), maxData(:,2), maxMod(:,2));
    fullMMTable.Properties.VariableNames = {'cellNum', 'iterNum', 'index', 'normFBPos', 'normSBPos', 'speedCor', 'timeDiff', 'width', 'duration', 'FBMaxData', 'FBMaxMod', 'SBMaxData', 'SBMaxMod'};
    clear tempTab maxData maxMod tempMM* tempT sbInd fbInd *bTime 
    
    allCellDat.fullMMTable = fullMMTable; 

end
    
totTab = vertcat(sbTable, mbTable, sgTable, mgTable, mmTable);
numRows = height(totTab); 
cellN = cellNum * ones(numRows, 1); 
iterN = iterNum * ones(numRows, 1); 


allCellDat.table= [table(cellN, iterN), totTab]; 
allCellDat.violTable = violTab; 


    
end




