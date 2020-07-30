function minMotSt = calcMinMotExtLinCompDiffWandV(mmStruct, sbStruct)


%function minMotSt = calcMinMotExtLinCompDiffWandV(pStruct)
%
% This function is designed specifically for minimal motion protocols that
% are diagonnaly corrected (with gratingTable and not gratingInds). It find
% the single bar presentations and generates the predicted linear sum by
% shifting them appropriately. 
%
% This function is a modification of calcMinMotExtLinCompDiff with the
% ability to take different widths and bar values into account. It also
% uses the singleBar input to add into the diagonal calculation when
% avaliable
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
%                   protocol. contains the following fields for each
%                   [firstBar, secondBar, timing] 
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


error('Bad version of the function! do not use!')


% base parameters

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation
barCombOrder = [0, 0; 1, 0; 0, 1]; % since the more relevant question is with preceding B (and oder is to remain consistent across cells
spCorrOrd = [1, 0]; % so that speed corrected will be the first position
% To be consistent across cells
timeDiffOrd = [0, 0.02, 0.04, 0.08, 0.16, 0.32]; 
widOrd = [1,2,4];

% checking inputs
relTab = mmStruct.gratingTable;

if any(relTab.FBStat == 1)
    relTab = relTab(relTab.FBStat == 0, :);
    warning('This funtion only calculates FBStat = 0, FBStat =1 stimuli not included in the resulting structure')
end

assert(all(ismember(['fAppear'; 'sAppear'], relTab.Properties.VariableNames)), ...
           'gratingTable is missing either fAppear or sAppear')

assert(ismember('FBStat', relTab.Properties.VariableNames), 'This function is designed for minMot protocols only')
assert(all(unique(relTab.FBStat) == 0), 'this function is designed for FBStat zero only')

ufPos = unique(relTab.FBPos);
usPos = unique(relTab.SBPos);
spCorr = unique(relTab.speedCor);

uMMPos = union(ufPos, usPos);

barW = unique(relTab.width);
uTD = unique(relTab.timeDiff);

barComb = relTab{:, {'FBVal', 'SBVal'}};
uniBarC = unique(barComb, 'rows');

singleBarSt = generateAlignedSingleBarStwMinMaxDiffWandV(sbStruct);
sbTab = sbStruct.gratingTable; 
sbPos = unique(sbTab.position);
sbDur = unique(sbTab.stimDur);
sbVal = unique(sbTab.value);
sbWid = unique(sbTab.width);

relMaxExt = singleBarSt(end, end, end, end).maxExtPosVal;
relPD = sign(singleBarSt(end, end, end, end).maxInhPosVal - relMaxExt);
    
genSBPos = (sbPos - relMaxExt)*relPD; % zeroed SBPos 

% assert(all(ufPos == usPos), 'FBPos must be equal to SBPos') Not sure this
% is necessary

% alignSt = alignProtocolDataByTable(pStruct, 'sAppear'); % not compatible
% with later functions
alignSt = alignProtocolDataByTable(mmStruct, 'fAppear');

% orginizing mean data

minMotSt = struct;


for cc=1:length(spCorr)
    
    relSPCInd = find(spCorrOrd == spCorr(cc));

    for ww=1:length(barW)
        
        widI = find(widOrd == barW(ww)); 

        for bc=1:size(uniBarC,1)

            relComb = uniBarC(bc, :); 
            
            bcI = find(all(barCombOrder == relComb, 2));

            FBVal = relComb(1);
            SBVal = relComb(2);

            for fp=1:length(ufPos)
                
                fpInd = find(uMMPos == ufPos(fp));
                zFBPos = genSBPos(sbPos == ufPos(fp));

                for sp=1:length(usPos)
                    
                    spInd = find(uMMPos == usPos(sp));
                    zSBPos = genSBPos(sbPos == usPos(sp));

                    for tt=1:length(uTD)
                        
                        tdI = find(timeDiffOrd == uTD(tt)); 

                        relInd = ismember(relTab{:, {'FBVal'; 'SBVal'; 'FBPos'; 'SBPos'; 'timeDiff'; 'width'; 'speedCor'}}, ...
                                          [FBVal, SBVal, ufPos(fp), usPos(sp), uTD(tt), barW(ww), spCorr(cc)], 'rows');

                        if unique(relInd) == 0
                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).empty = 1; 
                            continue
                        else

                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).data = alignSt(relInd);

                            relDat = minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).data.align.mean;
                            baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
                            zeroInd = find(relDat(:,1) > 0, 1, 'first');

                            baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
                            baseSubResp = relDat;
                            baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).subData.baseSub = baseSubResp;
                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).subData.baseline = baseVal;
                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).subData.zeroInd = zeroInd;
                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).subData.length = size(baseSubResp, 1);
                            minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).empty = 0;
%                             tempTab = minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).data.table; 
%                             minMotSt(relSPCInd, widI, bcI, fpInd, spInd, tdI).data.table = [tempTab, table(zFBPos, zSBPos)];

                        end

                    end

                end

            end

        end

    end
       
end

% a bit silly but didn't have a better solution
% fill in all the combination that didn't even appear in the above for
% loops

dSiz = size(minMotSt);

for ii=1:dSiz(1)
    for jj=1:dSiz(2)
        for kk=1:dSiz(3)
            for mm=1:dSiz(4)
                for nn=1:dSiz(5)
                    for tt=1:dSiz(6)
                        if isempty(minMotSt(ii,jj,kk,mm,nn,tt).empty)
                            minMotSt(ii,jj,kk,mm,nn,tt).empty = 1; 
                        end
                    end
                end
            end
        end
    end
end


% verify protocols have the same input coordinates
sbCen = sbStruct.inputParams.gridCenter;
mmCen = mmStruct.inputParams.gridCenter;
mmTimeDiffZeroDur = mmStruct.inputParams.stepDur; % duration when timeDiff is 0 

assert(all(sbCen == mmCen), 'center are not identical - add offset')

for cc=1:length(spCorr)
    
    relSPCInd = find(spCorrOrd == spCorr(cc));

    for ww=1:length(barW)
        
        widI = find(widOrd == barW(ww)); 
        sbWI = find(sbWid == barW(ww)); 
        
        for bc=1:size(uniBarC,1)

            relComb = uniBarC(bc, :); 
            bcI = find(all(barCombOrder == relComb, 2));

            FBVal = relComb(1);
            SBVal = relComb(2);
            
            sbFVI = find(sbVal == FBVal); 
            sbSVI = find(sbVal == SBVal); 

            for tt=1:length(uTD)
                
                relT = uTD(tt);
                tdI = find(timeDiffOrd == relT);
                
                if relT == 0 
                    relT = mmTimeDiffZeroDur;
                end
                sbTDI = find(sbDur == relT);

                for fp=1:length(ufPos)
                    
                    fpI = find(uMMPos == ufPos(fp));
                    
                    if fpI <= dSiz(5)
                        fbTemp = minMotSt(relSPCInd, widI, bcI, fpI, fpI, tdI);
                    else
                        fbTemp.empty = 1; 
                    end
                    
                    sbFPosI = find(sbPos == ufPos(fp));
                    
                    % if there is relevant singleBar data, combine with
                    % minMot
                   if length([sbFPosI, sbTDI, sbWI, sbFVI]) == 4 % all indices are viable (still doesn't mean it has data)
                       sbFbTemp = singleBarSt(sbFPosI, sbTDI, sbWI, sbFVI);  
                   else
                       sbFbTemp.empty = 1; 
                   end
                        
                   if sbFbTemp.empty && fbTemp.empty
                       continue
                   elseif sbFbTemp.empty && ~fbTemp.empty
                       combFBSt = fbTemp.subData;
                   elseif ~sbFbTemp.empty && fbTemp.empty
                       combFBSt = sbFbTemp.subData;
                   else % if both have data
                        % finding which one is longer and combining them
                        
                        numRepF = length(fbTemp.data.align.rep);
                        numRepSBF = length(sbFbTemp.data.align.rep);
                        wRF = numRepF/(numRepF+numRepSBF);
                        wRSBF = numRepSBF/(numRepF+numRepSBF);
                        
                        if sbFbTemp.subData.length > fbTemp.subData.length
                            
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
                        combFBSt.baseSub(:,2) = shiftVec * wShift + nonShiftedVec(:,2) * wNShift; 
                        
                   end     
                        

                    for sp=1:length(usPos)
                        
                        spI = find(uMMPos == usPos(sp));
                        
                        if spI <= dSiz(4)
                            sbTemp = minMotSt(relSPCInd, widI, bcI, spI, spI, tdI);
                        else
                            sbTemp.empty = 1; 
                        end
                        
                        sbSPosI = find(sbPos == usPos(sp));
                    
                        % if there is relevant singleBar data, combine with
                        % minMot
                       if length([sbSPosI, sbTDI, sbWI, sbSVI]) == 4 % all indices are viable (still doesn't mean it has data)
                           sbSbTemp = singleBarSt(sbSPosI, sbTDI, sbWI, sbSVI);  
                       else
                           sbSbTemp.empty = 1; 
                       end
                       
                       if sbSbTemp.empty && sbTemp.empty
                        continue
                       elseif sbSbTemp.empty && ~sbTemp.empty
                           combSBSt = sbTemp.subData;
                       elseif ~sbSbTemp.empty && sbTemp.empty
                           combSBSt = sbSbTemp.subData;
                       else % if both have data
                            % finding which one is longer and combining them
                            numRepS = length(sbTemp.data.align.rep);
                            numRepSBS = length(sbSbTemp.data.align.rep);
                            wRS = numRepS/(numRepS+numRepSBS);
                            wRSBS = numRepSBS/(numRepS+numRepSBS);
                        
                            if sbSbTemp.subData.length > sbTemp.subData.length

                                combSBSt = sbSbTemp.subData; 

                                relLen = sbSbTemp.subData.length;
                                newZero = sbSbTemp.subData.zeroInd; 
                                inputSt.zeroInd = sbTemp.subData.zeroInd;
                                inputSt.data = sbTemp.subData.baseSub(:,2);
                                shiftVec = padRespVecGen(inputSt, newZero, relLen); 
                                nonShiftedVec = sbSbTemp.subData.baseSub; % with time
                                wShift = wRS;
                                wNShift = wRSBS; 

                            else

                                combSBSt = sbTemp.subData; 

                                relLen = sbTemp.subData.length;
                                newZero = sbTemp.subData.zeroInd; 
                                inputSt.zeroInd = sbSbTemp.subData.zeroInd;
                                inputSt.data = sbSbTemp.subData.baseSub(:,2);
                                shiftVec = padRespVecGen(inputSt, newZero, relLen); 
                                nonShiftedVec = sbTemp.subData.baseSub; % with time

                                wShift = wRSBS;
                                wNShift = wRS; 

                            end

                            % weighted average of singleBar and minMot single
                            % response (weights are repeats)
                            combSBSt.baseSub(:,2) = shiftVec * wShift + nonShiftedVec(:,2) * wNShift; 

                       end     
                       

                        if minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).empty
                            linData = 0;
                            minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).data.table = table(relSPCInd, widI, bcI, fpI, spI, tdI, linData);
                            continue
                            
                        elseif fpI == spI
                            linData = 0;
                            tempTab = minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).data.table;
                            minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).data.table = [tempTab, table(relSPCInd, widI, bcI, fpI, spI, tdI, linData)];
                            continue
                            
                        else

                            % instead of estimating, this actually takes the mean
                            % position for the second bar appearance and find the
                            % appropriate time to copy the linSum response to
                            relSt = minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI); 
                            totLen = relSt.subData.length;
                            relZeroInd = relSt.subData.zeroInd;
                            sApFr = relSt.data.table.sAppear;
                            oriSApInd = relSt.data.align.meanPos(relSt.data.align.meanPos(:,2) == sApFr, 1);
                            
                            inputStF.zeroInd = combFBSt.zeroInd;
                            inputStF.data = combFBSt.baseSub(:,2);
                            shiftFVec = padRespVecGen(inputStF, relZeroInd, totLen); 
                            
                            inputStS.zeroInd = combSBSt.zeroInd;
                            inputStS.data = combSBSt.baseSub(:,2);
                            shiftSVec = padRespVecGen(inputStS, oriSApInd, totLen); 
                            
                            minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).linSum = [relSt.subData.baseSub(:,1), shiftFVec + shiftSVec];
                            minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).fbVec = shiftFVec;  % for easier QC
                            minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).sbVec = shiftSVec; 
                            tempTab = minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).data.table;
                            linData = 1; 
                            minMotSt(relSPCInd, widI, bcI, fpI, spI, tdI).data.table = [tempTab, table(relSPCInd, widI, bcI, fpI, spI, tdI, linData)];
                        end

                    end

                end

            end
            
        end
        
    end
    
end




                

end
       
       
       
       