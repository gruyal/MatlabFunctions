function minMotSt = calcMinMotExtLinCompDiffWandVwTable(mmStruct, sbStruct)

%function minMotSt = calcMinMotExtLinCompDiffWandV(pStruct)
%
% This function is designed specifically for minimal motion protocols that
% are diagonnaly corrected (with gratingTable and not gratingInds). It find
% the single bar presentations and generates the predicted linear sum by
% shifting them appropriately. 
%
% This function is a modification of calcMinMotExtLinCompDiff, however it 
% does not organize the with the data in a arrayed structue, but instead
% leaves it in the same order as in the aligned structure and adds a table
% to it for reference. This is done since the new round of experiments
% created a 6D structure which was unwielding 
%
% function is designed to work on combMinMotStruct, which is the structure
% of all the minMot protocols from a single cell combined. 
%
%
% INPUT 
%
% mmStruct -        minMot protocol with manual addition of bar
%                   presentation frames to gratingTable. 
%                   specifically fAppear and sAppear 
% sbStruct =        singleBar protocol from the same cell
%
% OUTPUT
% minMotSt -        structure that organizes the output from the minMot
%                   protocol. contains the following fields:
% .mmTable -         Table from mmStruct. Data is organized in minMotSt in
%                   the same order as this table. Table is also added
%                   several variable:
%                   normFB/SBPos - normalized position after PD and center
%                   have been calculated
%                   normPosDiff - difference between above 2 vars
%                   missBar - for combination where sum have not been
%                   computed mentions which bar was missing (0,1,2)
%                   sbFB/SBIndex - index from sbTable for the data that was used for the individual bar position 
%
% .sbTable -        table from sbStruct. included as a reference for QC (to
%                   see which stim were used for the calcultion of which
%                   combination)
% .mmResult -       structure with the following fields where indices are
%                   indentical to mmTable
%
%  .data -
%    .align -       Generate by alignProtocolDataByTable while alinging to
%                   fAppear. Contain the following fields
%       .rep -      all the individual repeats for the stim
%       .mean -     mean response
%       .meanPos-   mean position and timing matrix (timimng corresponds to
%                   the mean fieid
%    .table -       relevant table row for the above data
%
%  .subData -       baseline subtracted data. Contains the following fields
%     .baseline -   the calculated baseline that was subtracted
%     .zeroInd -    index to when the time stamp crosses zero
%     .length -     number of samples in data
%     .baseSub -    baseline subtracted data
%  .linSum -        if applicable timing and data vector constructed by a
%                   linear summation of the individual bar responses
%  .fb/sbVec -      shifted vectors from singleBar data (for QC)
%  .sbInd -         index of when the second bar appears and disappears (FB is zero)
%
% .combSingBarSt -  added this structure as QC to see how the singlebar
%                   responses are being combined. Has the idividual vectors
%                   their weights and their combination. 



% base parameters

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

% checking inputs
mmTab = mmStruct.gratingTable; 

if any(mmTab.FBStat == 1)
    warning('This funtion only calculates FBStat = 0, FBStat = 1 stimuli will be wrongfully calculated')
end

assert(all(ismember(['fAppear'; 'sAppear'], mmTab.Properties.VariableNames)), ...
           'gratingTable is missing either fAppear or sAppear')

assert(ismember('FBStat', mmTab.Properties.VariableNames), 'This function is designed for minMot protocols only')
% assert(all(unique(mmTab.FBStat) == 0), 'this function is designed for FBStat zero only')

singleBarSt = generateAlignedSingleBarStwMinMaxDiffWandV(sbStruct);
sbTab = sbStruct.gratingTable; 
sbPos = unique(sbTab.position);
sbDur = unique(sbTab.stimDur);
sbVal = unique(sbTab.value);
sbWid = unique(sbTab.width);

relMaxExt = singleBarSt(end, end, end, end).maxExtPosVal;
relPD = sign(singleBarSt(end, end, end, end).maxInhPosVal - relMaxExt);

% alignSt = alignProtocolDataByTable(pStruct, 'sAppear'); % not compatible
% with later functions
alignSt = alignProtocolDataByTable(mmStruct, 'fAppear');
allMMStim = length(alignSt);

% orginizing mean data

minMotSt = struct;

sbTab = sbStruct.gratingTable; 
minMotSt.sbTable = sbTab; 


for ii=1:allMMStim
    
    minMotSt.mmResult(ii).data = alignSt(ii); 
    relDat = minMotSt.mmResult(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    minMotSt.mmResult(ii).subData.baseSub = baseSubResp;
    minMotSt.mmResult(ii).subData.baseline = baseVal;
    minMotSt.mmResult(ii).subData.zeroInd = zeroInd;
    minMotSt.mmResult(ii).subData.length = size(baseSubResp, 1);
    
end
    
    


% verify protocols have the same input coordinates
sbCen = sbStruct.inputParams.gridCenter;
mmCen = mmStruct.inputParams.gridCenter;
mmTimeDiffZeroDur = mmStruct.inputParams.stepDur; % duration when timeDiff is 0 

% assert(all(sbCen == mmCen), 'centers are not identical - add offset')

 if any(sbCen ~= mmCen) 
     offset = input('centers are not identical - add offset \n');
     assert(round(offset) == offset, 'offset should be an integer')
 else
     offset = 0;
 end

missBar = zeros(allMMStim, 1); 
sbFBIndex = zeros(allMMStim, 1); 
sbSBIndex = zeros(allMMStim, 1); 

count=0;
jCount = 0;
justSingBarSt = struct; 
combSingBarSt = struct; 

for ii=1:allMMStim
    
    relTab = mmTab(ii, :);
    relW = relTab.width; 
    relFBV = relTab.FBVal; 
    relSBV = relTab.SBVal; 
    relFBP = relTab.FBPos; 
    relSBP = relTab.SBPos;
    relTD = relTab.timeDiff; 
    relSP = relTab.speedCor; 
    
    if relFBP == relSBP
        continue
    end
    
    % since in later protocols (after March 2018) diagonal presented flicker
    % FB SB if values were not the same. Before presented only second bar
    if relFBV == relSBV 
        diagFBInd = ismember(mmTab{:, {'FBVal'; 'SBVal'; 'FBPos'; 'SBPos'; 'timeDiff'; 'width'; 'speedCor'}}, ...
                                          [relFBV,relSBV, relFBP, relFBP, relTD, relW, relSP], 'rows'); 
        diagSBInd = ismember(mmTab{:, {'FBVal'; 'SBVal'; 'FBPos'; 'SBPos'; 'timeDiff'; 'width'; 'speedCor'}}, ...
                                          [relFBV,relSBV, relSBP, relSBP, relTD, relW, relSP], 'rows');
        
        fbTemp = minMotSt.mmResult(diagFBInd); 
        sbTemp = minMotSt.mmResult(diagSBInd); 
        
    else
        fbTemp = [];
        sbTemp = [];
        
    end
    
    % when 2 bars are presented sim, there is an additional paramter to
    % determine duration
    if relTD == 0
        relTD = mmTimeDiffZeroDur;
    end
    
    relSbFPos = sbPos == relFBP + offset; 
    relSbSPos = sbPos == relSBP + offset; 
    relSbTD = sbDur == relTD; 
    relSbW = sbWid == relW;
    relSbFBV = sbVal == relFBV; 
    relSbSBV = sbVal == relSBV; 
    
    sbFbTemp = singleBarSt(relSbFPos, relSbTD, relSbW, relSbFBV); 
    sbSbTemp = singleBarSt(relSbSPos, relSbTD, relSbW, relSbSBV); 
    
    
    %combining FB data from minMot diagonal and singleBar 
    if (isempty(sbFbTemp) || sbFbTemp.empty) && isempty(fbTemp)
       missBar(ii) = 1; 
       continue
    elseif (isempty(sbFbTemp) || sbFbTemp.empty) && ~isempty(fbTemp)
        combFBSt = fbTemp.subData;
    elseif ~sbFbTemp.empty && isempty(fbTemp)
        
        combFBSt = sbFbTemp.subData;
        
        % add a new enetry only if it isn't there already 
        if ~ismember(sbFbTemp.data.table.index, sbFBIndex)
            jCount=jCount+1;
            justSingBarSt(jCount).sbTabIndex = sbFbTemp.data.table.index; 
            justSingBarSt(jCount).singleBarVec = combFBSt.baseSub(:,2);
        end
        
        sbFBIndex(ii) = sbFbTemp.data.table.index;
        
    else % if both have data
        % finding which one is longer and combining them

        numRepF = length(fbTemp.data.align.rep) - length(fbTemp.data.exclude);
        numRepSBF = length(sbFbTemp.data.align.rep) - length(sbFbTemp.data.exclude);
        wRF = numRepF/(numRepF+numRepSBF);
        wRSBF = numRepSBF/(numRepF+numRepSBF);
        
        if sbFbTemp.subData.length > fbTemp.subData.length
            
            sbLongFlag = 1;
            combFBSt = sbFbTemp.subData; 

            relLen = sbFbTemp.subData.length;
            newZero = sbFbTemp.subData.zeroInd; 
            inputSt.zeroInd = fbTemp.subData.zeroInd;
            inputSt.data = fbTemp.subData.baseSub(:,2);
            shiftVec = padRespVecGen(inputSt, newZero, relLen); 
            nonShiftedVec = sbFbTemp.subData.baseSub; % with time
            wShift = wRF;
            wNShift = wRSBF; 
        
        else
            sbLongFlag = 2;
            combFBSt = fbTemp.subData;

            relLen = fbTemp.subData.length;
            newZero = fbTemp.subData.zeroInd; 
            inputSt.zeroInd = sbFbTemp.subData.zeroInd;
            inputSt.data = sbFbTemp.subData.baseSub(:,2);
            shiftVec = padRespVecGen(inputSt, newZero, relLen); 
            nonShiftedVec = fbTemp.subData.baseSub; % with time

            wShift = wRSBF;
            wNShift = wRF; 

        end
        % weighted average of singleBar and minMot single
        % response (weights are repeats)
        combFBSt.baseSub(:,2) = nansum([shiftVec * wShift, nonShiftedVec(:,2) * wNShift], 2); 
        tempVecs = [nonShiftedVec(:,2), shiftVec];
        tempWs = [wNShift, wShift];
        
        % add a new enetry only if it isn't there already 
        if ~ismember(sbFbTemp.data.table.index, sbFBIndex)
            count=count+1;
            combSingBarSt(count).sbTabIndex = sbFbTemp.data.table.index; 
            combSingBarSt(count).singleBarVec = tempVecs(:,sbLongFlag);
            combSingBarSt(count).mmDiagVec = tempVecs(:,3-sbLongFlag);
            combSingBarSt(count).combined = combFBSt.baseSub(:,2);
            combSingBarSt(count).singleBarW = tempWs(sbLongFlag); 
            combSingBarSt(count).mmDiagW = tempWs(3-sbLongFlag);
        end
        
        sbFBIndex(ii) = sbFbTemp.data.table.index; 
    end
    
    %combining SB data from minMot diagonal and singleBar 
    if (isempty(sbSbTemp) || sbSbTemp.empty) && isempty(sbTemp)
       missBar(ii) = 1; 
       continue
    elseif (isempty(sbSbTemp) || sbSbTemp.empty) && ~isempty(sbTemp)
        combSBSt = sbTemp.subData;
    elseif ~sbSbTemp.empty && isempty(sbTemp)
        sbSBIndex(ii) = sbSbTemp.data.table.index; 
        combSBSt = sbSbTemp.subData;
    else % if both have data
        % finding which one is longer and combining them
        
        sbSBIndex(ii) = sbSbTemp.data.table.index; 
        numRepF = length(sbTemp.data.align.rep) - length(sbTemp.data.exclude);
        numRepSBF = length(sbSbTemp.data.align.rep) - length(sbSbTemp.data.exclude);
        wRF = numRepF/(numRepF+numRepSBF);
        wRSBF = numRepSBF/(numRepF+numRepSBF);
        
        if sbSbTemp.subData.length > sbTemp.subData.length
            
            combSBSt = sbSbTemp.subData; 

            relLen = sbSbTemp.subData.length;
            newZero = sbSbTemp.subData.zeroInd; 
            inputSt.zeroInd = sbTemp.subData.zeroInd;
            inputSt.data = sbTemp.subData.baseSub(:,2);
            shiftVec = padRespVecGen(inputSt, newZero, relLen); 
            nonShiftedVec = sbSbTemp.subData.baseSub; % with time
            wShift = wRF;
            wNShift = wRSBF; 
        
        else
            
            combSBSt = sbTemp.subData;

            relLen = sbTemp.subData.length;
            newZero = sbTemp.subData.zeroInd; 
            inputSt.zeroInd = sbSbTemp.subData.zeroInd;
            inputSt.data = sbSbTemp.subData.baseSub(:,2);
            shiftVec = padRespVecGen(inputSt, newZero, relLen); 
            nonShiftedVec = sbTemp.subData.baseSub; % with time

            wShift = wRSBF;
            wNShift = wRF; 

        end
        % weighted average of singleBar and minMot single
        % response (weights are repeats)
        combSBSt.baseSub(:,2) = nansum([shiftVec * wShift, nonShiftedVec(:,2) * wNShift], 2); 
    
    end
    
    
    % Creating linear comparison
    relSt = minMotSt.mmResult(ii); 
    totLen = relSt.subData.length;
    relZeroInd = relSt.subData.zeroInd;
    sApFr = relSt.data.table.sAppear;
    sDisFr = relSt.data.table.sDisappear;
    oriSApInd = relSt.data.align.meanPos(relSt.data.align.meanPos(:,2) == sApFr, 1);
    oriSDisInd = relSt.data.align.meanPos(relSt.data.align.meanPos(:,2) == sDisFr, 1);

    inputStF.zeroInd = combFBSt.zeroInd;
    inputStF.data = combFBSt.baseSub(:,2);
    shiftFVec = padRespVecGen(inputStF, relZeroInd, totLen); 

    inputStS.zeroInd = combSBSt.zeroInd;
    inputStS.data = combSBSt.baseSub(:,2);
    shiftSVec = padRespVecGen(inputStS, oriSApInd, totLen); 

    minMotSt.mmResult(ii).linSum = [relSt.subData.baseSub(:,1), nansum([shiftFVec, shiftSVec], 2)];
    minMotSt.mmResult(ii).fbVec = shiftFVec;  % for easier QC
    minMotSt.mmResult(ii).sbVec = shiftSVec; 
    minMotSt.mmResult(ii).sbInd = [oriSApInd, oriSDisInd]; 
    
end


minMotSt.mmTable = mmTab;

if relPD == 1
    normFBPos = (minMotSt.mmTable.FBPos - relMaxExt) * relPD; 
    normSBPos = (minMotSt.mmTable.SBPos - relMaxExt) * relPD;
else
    normFBPos = (minMotSt.mmTable.FBPos - relMaxExt - minMotSt.mmTable.width + 1) * relPD; 
    normSBPos = (minMotSt.mmTable.SBPos - relMaxExt - minMotSt.mmTable.width + 1) * relPD;
end
    
normPosDiff = normSBPos - normFBPos; 
minMotSt.mmTable = [minMotSt.mmTable, table(normFBPos, normSBPos, normPosDiff, missBar, sbFBIndex, sbSBIndex)];

minMotSt.combSingBarSt = combSingBarSt;
minMotSt.justSingBarSt = justSingBarSt;



                

end
       
       
       
       